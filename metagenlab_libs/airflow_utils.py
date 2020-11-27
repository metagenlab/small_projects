

import os
import yaml
from datetime import datetime
from metagenlab_libs import gendb_utils


def clean_species(species_string):
    import re 
    # case when fastq file could not be matched to sample table
    if species_string is None:
        return "Unspecified species"
    if len(species_string) == 0:
        return "Unspecified species"
    species = re.sub(",.*", "",species_string)
    species = " ".join(species.split(" ")[0:2])
    return species 



def make_run_dir(execution_folder, 
                 analysis_id):
    import shutil
    import errno
    
    folder_name = datetime.now().strftime("%Y_%m_%d-%H%M")
    
    run_execution_folder = os.path.join(execution_folder, folder_name)
    
    if eval(analysis_id):
        print("adding execution folder for analysis", analysis_id)
        print(type(analysis_id))
        gendb_utils.add_analysis_metadata(analysis_id, "airflow_execution_folder", run_execution_folder)
    
    # remove existing dir
    shutil.rmtree(run_execution_folder, ignore_errors=True)
    # make dir
    os.makedirs(run_execution_folder)
    
    return folder_name

def write_sample_file(gen_db,
                      fastq_list,
                      analysis_id, 
                      execution_folder):
        
    if not isinstance(fastq_list, list):
        fastq_list = fastq_list.split(",")
    
    # id,fastq_prefix,R1,R2,species_name
    fastq_df = gen_db.get_fastq_metadata(fastq_list)
    
    run_execution_folder = os.path.join(execution_folder, analysis_id)
    
    header = ["SampleName",
              "ScientificName",
              "R1",
              "R2", 
              "fastq_id"]
    
    with open(os.path.join(run_execution_folder, f'{analysis_id}.tsv'), 'w') as f:
        f.write("\t".join(header) + '\n')
        for n, row in fastq_df.iterrows():
            
            # deal with multiple fastq with same name
            R1 = row["R1"]
            R2 = row["R2"]
            fastq_id = row["fastq_id"]
            species = row["species_name"]
            sample_name = f'{row["fastq_prefix"]}_{fastq_id}'
            f.write(f"{sample_name}\t{species}\t{R1}\t{R2}\t{fastq_id}\n")
            
            
            
def write_snakemake_config_file(analysis_id,
                                fastq_list,
                                execution_folder,
                                snakemake_config,
                                gen_db,
                                reference_list=False,
                                scientific_name=False,
                                check_single_species=False):
    
    run_execution_folder = os.path.join(execution_folder, analysis_id)
    
    if check_single_species and not scientific_name:
        species_list = list(set(gen_db.get_fastq_id2species(fastq_list.split(",")).values()))
    if check_single_species:
        if len(species_list) > 1:
            raise IOError("More than one different species in the dataset: %s" % ','.join(species_list))
        else:
            scientific_name = species_list[0]
    
    # if references, prepare list
    if reference_list:
        reference_fastq_list = reference_list.split(",")
        fastq_df = gen_db.get_fastq_metadata(reference_fastq_list)
        ref_list = []
        for n, row in fastq_df.iterrows():
            fastq_id = row["fastq_id"]
            sample_name = f'{row["fastq_prefix"]}_{fastq_id}'
            ref_list.append(sample_name)
    
    with open(os.path.join(run_execution_folder, f'{analysis_id}.config'), 'w') as f:
        # update sample table name
        snakemake_config["local_samples"] = f'{analysis_id}.tsv'
        if reference_list:
            snakemake_config["reference"] = f'{",".join(ref_list)}'
        if check_single_species:
            snakemake_config["species"] = f'{scientific_name}'

        documents = yaml.dump(snakemake_config, f)