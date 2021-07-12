

import os
import yaml
from datetime import datetime
from metagenlab_libs import gendb_utils
from airflow.hooks.base_hook import BaseHook
from airflow.providers.slack.operators.slack_webhook import SlackWebhookOperator
import string 
import random
import gzip
import shutil
  
SLACK_CONN_ID = 'slack'

def compress(input, output):
    with open(input, 'rb') as f_in:
        if not output.endswith(".gz"):
            output += '.gz'
        with gzip.open(output, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
def list_dir(root_dir):
    file_set = set()
    for dir_, _, files in os.walk(root_dir):
        for file_name in files:
            rel_dir = os.path.relpath(dir_, root_dir)
            rel_file = os.path.join(rel_dir, file_name)
            file_set.add(rel_file)
    return file_set

def copy_and_compress(source, target, compress_ext):
    '''
    Input:
    - source directory and target directory 
    - source filename and target filename
    (no mix of file and directory)
    '''
    if os.path.isdir(source):
        # list of relative path
        complete_file_lst = list_dir(source)
        for filename in complete_file_lst:
            # if directory does not exist, create it
            source_abs_path = os.path.join(source, filename)
            target_abs_path = os.path.join(target, filename)
            if not os.path.exists(os.path.dirname(target_abs_path)):
                os.makedirs(os.path.dirname(target_abs_path))
            # check extension
            extension = filename.split(".")[-1]
            if extension in compress_ext:
                compress(source_abs_path, target_abs_path)
            else:
                shutil.copy(source_abs_path, target_abs_path)
    else:
        if not os.path.exists(os.path.dirname(target)):
            os.makedirs(os.path.dirname(target))
        extension = source.split(".")[-1]
        if extension in compress_ext:
            compress(source, target)
        else:
            shutil.copy(source, target)


def clean_species(species_strings):s
    import re 
    # case when fastq file could not be matched to sample table
    if species_string is None:
        return "Unspecified species"
    if len(species_string) == 0:
        return "Unspecified species"
    species = re.sub(",.*", "",species_string)
    species = " ".join(species.split(" ")[0:2])
    return species 


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def make_run_dir(execution_folder, 
                 analysis_id=False):
    import shutil
    import errno
    
    folder_name = datetime.now().strftime("%Y_%m_%d-%H%M_" + id_generator(6))
    
    run_execution_folder = os.path.join(execution_folder, folder_name)
    
    if analysis_id:
        # update status
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)
        # save execution folder
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_folder", run_execution_folder)
    
    # remove existing dir
    shutil.rmtree(run_execution_folder, ignore_errors=True)
    # make dir
    os.makedirs(run_execution_folder)
    
    return folder_name


def backup_output_files_samples(metadata_name2path_template, 
                                fastq_list,
                                analysis_name,
                                analysis_id,
                                backup_folder):
    
    GEN_DB = gendb_utils.DB()
    
    fastq_df = GEN_DB.get_fastq_and_sample_data(fastq_list)   
    
    qc_data = [] 
    for metadata_name in metadata_name2path_template:
        path_template = metadata_name2path_template[metadata_name]
        for n, sample in fastq_df.iterrows():

            sample_name = f'{sample["sample_name"]}_{sample["fastq_id"]}'
            
            # assume structure: {workflow}/{analysis_name}/{filepath}
            backup_path_format_relative = path_template.format(analysis_name=analysis_name,sample=sample_name)
            backup_path_format_absolute = os.path.join(backup_folder, '/'.join(backup_path_format_relative.split("/")[1:]))
            
            print("backup_path_format", backup_path_format_relative)
            if not os.path.exists(backup_path_format_absolute):
                print(f"WARNING: {backup_path_format_absolute} does not exit, skipping" )
                continue

            qc_data.append({"fastq_id": sample["fastq_id"],
                            "metrics_name": metadata_name,
                            "metrics_value": backup_path_format_relative,
                            "pipeline_version": ""})
    for qc in qc_data:
        GEN_DB.add_fastq_metadata(fastq_id=qc["fastq_id"],
                                  term_name=qc["metrics_name"],
                                  value=qc["metrics_value"],
                                  analysis_id=analysis_id)
    


def backup_output_files_analysis(metadata_name2path_template, 
                                 analysis_name,
                                 analysis_id,
                                 backup_folder,
                                 save_backup_folder=True):
    
    print("backup_folder", backup_folder)

    for metadata in metadata_name2path_template:
        
        # assume structure: {workflow}/{analysis_name}/{filepath}
        path_template_format_relative = metadata_name2path_template[metadata].format(analysis_name=analysis_name)
        path_template_format_absolute = os.path.join(backup_folder, '/'.join(path_template_format_relative.split("/")[1:]))
        
        print("backup_path_format", path_template_format_absolute)
        if not os.path.exists(path_template_format_absolute):
            print(f"WARNING: {path_template_format_absolute} does not exit, skipping" )
            continue
        
        gendb_utils.add_analysis_metadata(analysis_id, 
                                          metadata, 
                                          path_template_format_relative, 
                                          update=True)
    if save_backup_folder:
        gendb_utils.add_analysis_metadata(analysis_id, 
                                          "backup_folder", 
                                          backup_folder.split("/")[-2] + f'/{analysis_name}', 
                                          update=True)

def write_sample_file(gen_db_path,
                      fastq_list,
                      analysis_name, 
                      execution_folder,
                      analysis_id):
    
    GEN_DB = gendb_utils.DB()
    
    # update status
    gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)
    
    if not isinstance(fastq_list, list):
        fastq_list = fastq_list.split(",")
    
    # id,fastq_prefix,R1,R2,species_name
    fastq_df = GEN_DB.get_fastq_and_sample_data(fastq_list)
    
    run_execution_folder = os.path.join(execution_folder, analysis_name)
    
    header = ["SampleName",
              "ScientificName",
              "R1",
              "R2", 
              "fastq_id"]
    
    with open(os.path.join(run_execution_folder, f'{analysis_name}.tsv'), 'w') as f:
        f.write("\t".join(header) + '\n')
        for n, row in fastq_df.iterrows():
            
            # deal with multiple fastq with same name
            R1 = row["R1"]
            R2 = row["R2"]
            fastq_id = row["fastq_id"]
            species = row["species_name"]
            sample_name = f'{row["sample_name"]}_{fastq_id}'
            f.write(f"{sample_name}\t{species}\t{R1}\t{R2}\t{fastq_id}\n")
            
            
            
def write_snakemake_config_file(analysis_name,
                                fastq_list,
                                execution_folder,
                                snakemake_config,
                                gen_db_path,
                                analysis_id,
                                reference_list=False,
                                check_single_species=False,
                                reference_docx=False,
                                additional_args=False):

    GEN_DB = gendb_utils.DB()
    
    # update status
    gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)

    run_execution_folder = os.path.join(execution_folder, analysis_name)
    
    species_list = list(set(GEN_DB.get_fastq_id2species(fastq_list.split(",")).values()))
    print("species_list", species_list)
    if check_single_species:
        if len(species_list) > 1:
            raise IOError("More than one different species in the dataset: %s" % ','.join(species_list))

    # if only one species, set scientific_name
    # otherwise Mixed
    if len(species_list) == 1:
        scientific_name = species_list[0]
    else:
        scientific_name = 'Mixed'
    
    print("reference list:", reference_list)
    # if references, prepare list
    if reference_list:
        reference_list = reference_list.split(",")
        fastq_df = GEN_DB.get_fastq_and_sample_data(reference_list)
        # check if external ref
        ref_list = [ref for ref in reference_list if ref not in fastq_df["fastq_id"].to_list()]
        if len(ref_list) != 0:
            print(f"WARNING: extrenal reference genome -- {ref_list[0]} ")
        ref_list += [f'{row["sample_name"]}_{row["fastq_id"]}' for n, row in fastq_df.iterrows()]
        if 'cgMLST' in reference_list:
            ref_list.append("cgMLST")
        
    print("ref list", ref_list)
    print("additional_args", additional_args)
    with open(os.path.join(run_execution_folder, f'{analysis_name}.config'), 'w') as f:
        # update sample table name
        snakemake_config["local_samples"] = f'{analysis_name}.tsv'
        if reference_list:
            snakemake_config["reference"] = f'{",".join(ref_list)}'
        snakemake_config["species"] = f'{scientific_name}'
        if reference_docx:
            snakemake_config["reference_docx"] = f'{reference_docx}'
        if additional_args:
            for arg in additional_args:
                snakemake_config[arg] = additional_args[arg]

        documents = yaml.dump(snakemake_config, f)
        

def backup(execution_folder, 
           backup_folder,
           file_or_folder_list,
           analysis_id=False,
           analysis_name=False,
           fastq_list=False,
           output_selection=False,
           config=False,
           workflow_name=False,
           compress_ext=["fna", "faa", "gbk", "gbff", "vcf", "tsv", "csv", "gff"]):
    
    '''
    Analysis name: folder within execution_folder (generally execution date)
    Analysis_metadata: dictionnary of status to add to LIMS: 
        {"analysis_id" : <id>,
         "value": <value>}
    '''
    
    import shutil
    import glob
    if analysis_id:
        # update status
        print("analysis_id", analysis_id)
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "running", update=True)

    # copy files and folders to backup directory
    for output in file_or_folder_list:
        
        if isinstance(output, list):
            # copy files to specified target directory
            target_dir = os.path.join(backup_folder,analysis_name, output[1])
            file_list = glob.glob(os.path.join(execution_folder, analysis_name, output[0]))
            print("file_list", file_list)

            for one_file in file_list:
                print("original", one_file)
                target_abs_path = os.path.join(target_dir, one_file)
                print("target", target_abs_path)
                copy_and_compress(original, target_abs_path, compress_ext)
                
        elif isinstance(output, dict):
            import re
            GEN_DB = gendb_utils.DB()
            # {glob: samples/*/mapping/bwa/*_assembled_genome.bam, regex: .*/samples/(.*)/mapping/bwa/(.*)_assembled_genome.bam, vars: {1: 'sample', 2: 'reference'}, target: "mapping/{sample}-vs-{reference}.bam", term: bam_file}
            file_list = glob.glob(os.path.join(execution_folder, analysis_name, output["glob"]))
            print("glob file list:", file_list)
            for one_file in file_list:
                s = re.search(output["regex"], one_file)
                term2value = {output["vars"][index]:s.group(index) for index in output["vars"]}
                term2value.update({'analysis_name': analysis_name})
                target_format = output["target"].format_map(term2value)
                target_path_full = os.path.join(backup_folder, '/'.join(target_format.split("/")[1:]))
                # copy file to target location
                if not os.path.exists(os.path.dirname(target_path_full)):
                    os.makedirs(os.path.dirname(target_path_full))
                print("cp:", one_file, target_path_full)
                copy_and_compress(one_file, target_path_full, compress_ext)
                # save path in db
                if "term" in output:
                    fastq_id = term2value["sample"].split("_")[-1]
                    GEN_DB.add_fastq_metadata(fastq_id=fastq_id,
                                            term_name=output["term"],
                                            value=target_format,
                                            analysis_id=analysis_id)
        else:
            # copy identical path
            output = output.format(analysis_name=analysis_name)
            original = os.path.join(execution_folder, analysis_name, output)
            target = os.path.join(backup_folder, analysis_name, output)
            print("original", original)
            print("target", target)
            # copy and compress what can be compressed
            copy_and_compress(original, target, compress_ext)

    # save file paths into database
    # can be either nested dictionnaries or a single dictionnary
    if analysis_name:
        if output_selection:
            output_selection = output_selection.split(",")
            metadata_lst = [value for key, value in config["WORKFLOW"][workflow_name]["PIPELINE_OUTPUT"]["ANALYSIS"].items() if key in output_selection]
            analysis_metadata_name2template = {k: v for d in metadata_lst for k, v in d.items()}
            
            if fastq_list:
                metadata_lst = [value for key, value in config["WORKFLOW"][workflow_name]["PIPELINE_OUTPUT"]["INDIVIDUAL_SAMPLES"].items() if key in output_selection]
                sample_metadata_name2template = {k: v for d in metadata_lst for k, v in d.items()}
        else:
            analysis_metadata_name2template = config["WORKFLOW"][workflow_name]["PIPELINE_OUTPUT"]["ANALYSIS"]
            if fastq_list:
                sample_metadata_name2template = config["WORKFLOW"][workflow_name]["PIPELINE_OUTPUT"]["INDIVIDUAL_SAMPLES"]
        print("backup path analysis", analysis_metadata_name2template)
        backup_output_files_analysis(analysis_metadata_name2template, 
                                     analysis_name,
                                     analysis_id,
                                     backup_folder)
        if fastq_list:
            fastq_list = fastq_list.split(",")
            backup_output_files_samples(sample_metadata_name2template, 
                                        fastq_list,
                                        analysis_name,
                                        analysis_id,
                                        backup_folder)



def task_fail_slack_alert(context):

    dag_run = context.get('dag_run')
    analysis_id = dag_run.conf.get('analysis_id')
    
    
    if analysis_id:
        print("Fail, updating analysis status")
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "failed", update=True)
        gendb_utils.update_analysis_status(analysis_id, "failed")
    else:
        print("NO ANALYSE NUMBER")

    print("task_fail_slack_alert")
    slack_webhook_token = BaseHook.get_connection(SLACK_CONN_ID).password

    slack_msg = """
            :red_circle: Task Failed. 
            *Task*: {task}  
            *Dag*: {dag} 
            *Execution Time*: {exec_date}  
            *Log Url*: {log_url} 
            *Analysis*: {analysis_id}
            """.format(
            task=context.get('task_instance').task_id,
            dag=context.get('task_instance').dag_id,
            exec_date=context.get('execution_date'),
            log_url=context.get('task_instance').log_url,
            analysis_id=analysis_id,
        )
          
    print("message:\n", slack_msg)
            
    failed_alert = SlackWebhookOperator(
        task_id='slack_test',
        http_conn_id='slack',
        webhook_token=slack_webhook_token,
        message=slack_msg,
        username='airflow')

    return failed_alert.execute(context=context)


def task_success_slack_alert(context):
    
    '''
    Slack message + update analysis status
    '''

    dag_run = context.get('dag_run')
    analysis_id = dag_run.conf.get('analysis_id')

    if analysis_id:
        print("Success, updating db")
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_status", "success", update=True)
        gendb_utils.update_analysis_status(analysis_id, "success")

    print("slack message")
    slack_webhook_token = BaseHook.get_connection(SLACK_CONN_ID).password

    slack_msg = """
            :heavy_check_mark: Dag Success. 
            *Dag*: {dag} 
            *Execution Time*: {exec_date}  
            *Log Url*: {log_url}
            *Analysis ID*: {analysis_id}
            """.format(
            dag=context.get('task_instance').dag_id,
            ti=context.get('task_instance'),
            exec_date=context.get('execution_date'),
            log_url=context.get('task_instance').log_url,
            analysis_id=analysis_id
        )
    
    success_alert = SlackWebhookOperator(
        task_id='slack_test',
        http_conn_id='slack',
        webhook_token=slack_webhook_token,
        message=slack_msg,
        username='airflow')

    return success_alert.execute(context=context)

