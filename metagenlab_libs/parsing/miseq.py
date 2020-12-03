import os
import re
import datetime

def parse_sample_sheet(sheet_path):
    '''
    [Header]
    IEMFileVersion,4
    Investigator Name,Sebastien Aeby
    Experiment Name,20200213_MTB
    Date,13.02.2020
    Workflow,GenerateFASTQ
    Application,FASTQ Only
    Assay,Nextera XT
    Description,
    Chemistry,Amplicon

    [Reads]
    151
    151

    [Settings]
    ReverseComplement,0
    Adapter,CTGTCTCTTATACACATCT

    [Data]
    Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
    193253,1052,,,N702,CGTACTAG,S504,AGAGTAGA,,
    191060,1053,,,N703,AGGCAGAA,S504,AGAGTAGA,,
    1801092178,1054,,,N705,GGACTCCT,S504,AGAGTAGA,,
    1907164214,1058,,,N710,CGAGGCTG,S504,AGAGTAGA,,
    20190718_positive_control,T+,,,N711,AAGAGGCA,S508,CTAAGCCT,,
    '''
    
    run2data = {}

    
    with open(sheet_path, 'r') as f:
        table_section = False
        reads_section = False
        sample_table = []
        read_len = []
        for row in f:
            if row.startswith("Date"):
                run_date = row.rstrip().split(",")[1]
                run_date_split = re.split("\.|\/", run_date)
                run_date = f'{int(run_date_split[2])}-{int(run_date_split[1]):02d}-{int(run_date_split[0]):02d}'
            if row.startswith("Experiment Name"):
                run_name = row.rstrip().split(",")[1]   
            if row.startswith("Assay"):
                assay = row.rstrip().split(",")[1]
            if row.strip() == '[Reads]':
                reads_section = True
                continue
            if row.strip() == '[Data]':
                table_section = True
                continue
            # if empty line or start of a new section
            if row.strip() == '' or row.strip()[0] == '[':
                table_section = False
                reads_section = False
            if table_section:
                sample_data = row.rstrip().split(",")
                sample_table.append(sample_data)
            if reads_section:
                read_len.append(row.rstrip())
    
    run2data["run_date"] = run_date
    run2data["run_name"] = run_name
    run2data["assay"] = assay
    run2data["reads"] = read_len
    run2data["path"] = os.path.dirname(sheet_path)
    
    return run2data


def parse_runinfo(runinfo_path):
    '''
<?xml version="1.0"?>
<RunParameters xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EnableCloud>true</EnableCloud>
  <RunParametersVersion>MiSeq_1_1</RunParametersVersion>
  <CopyManifests>true</CopyManifests>
  <FlowcellRFIDTag>
    <SerialNumber>000000000-BY7WY</SerialNumber>
    <PartNumber>15028382</PartNumber>
    <ExpirationDate>2019-04-26T00:00:00</ExpirationDate>
  </FlowcellRFIDTag>
  <PR2BottleRFIDTag>
    <SerialNumber>MS6813532-00PR2</SerialNumber>
    <PartNumber>15041807</PartNumber>
    <ExpirationDate>2019-05-08T00:00:00</ExpirationDate>
  </PR2BottleRFIDTag>
  <ReagentKitRFIDTag>
    <SerialNumber>MS6925444-500V2</SerialNumber>
    <PartNumber>15033573</PartNumber>
    <ExpirationDate>2019-04-03T00:00:00</ExpirationDate>
  </ReagentKitRFIDTag>
  <Resumable>true</Resumable>
  <ManifestFiles />
  <AfterRunWashMethod>Post-Run Wash</AfterRunWashMethod>
  <Setup>
    <SupportMultipleSurfacesInUI>true</SupportMultipleSurfacesInUI>
    <ApplicationVersion>2.6.2.1</ApplicationVersion>
    <NumTilesPerSwath>14</NumTilesPerSwath>
    <NumSwaths>1</NumSwaths>
    <NumLanes>1</NumLanes>
    <ApplicationName>MiSeq Control Software</ApplicationName>
  </Setup>
  <RunID>180823_M03935_0099_000000000-BY7WY</RunID>
  <ScannerID>M03935</ScannerID>
  <RunNumber>99</RunNumber>
  <FPGAVersion>9.5.12</FPGAVersion>
  <MCSVersion>2.6.2.1</MCSVersion>
  <RTAVersion>1.18.54</RTAVersion>
  <Barcode>000000000-BY7WY</Barcode>
  <PR2BottleBarcode>MS6813532-00PR2</PR2BottleBarcode>
  <ReagentKitPartNumberEntered>15033573</ReagentKitPartNumberEntered>
  <ReagentKitVersion>Version2</ReagentKitVersion>
  <ReagentKitBarcode>MS6925444-500V2</ReagentKitBarcode>
  <PreviousPR2BottleBarcode />
  <PreviousReagentKitBarcode />
  <ExperimentName>20180823_divers</ExperimentName>
  <Chemistry>Amplicon</Chemistry>
  <Username>sbsuser</Username>
  <Workflow>
    <Analysis>GenerateFASTQ</Analysis>
  </Workflow>
  <EnableAnalysis>false</EnableAnalysis>
  <Reads>
    <RunInfoRead Number="1" NumCycles="251" IsIndexedRead="N" />
    <RunInfoRead Number="2" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="3" NumCycles="8" IsIndexedRead="Y" />
    <RunInfoRead Number="4" NumCycles="251" IsIndexedRead="N" />
  </Reads>
  <TempFolder>D:\Illumina\MiSeqTemp\180823_M03935_0099_000000000-BY7WY</TempFolder>
  <AnalysisFolder>D:\Illumina\MiSeqAnalysis\180823_M03935_0099_000000000-BY7WY</AnalysisFolder>
  <RunStartDate>180823</RunStartDate>
  <MostRecentWashType>PostRun</MostRecentWashType>
  <RecipeFolder>D:\Illumina\MiSeq Control Software\CustomRecipe</RecipeFolder>
  <ILMNOnlyRecipeFolder>C:\Illumina\MiSeq Control Software\Recipe</ILMNOnlyRecipeFolder>
  <SampleSheetName>MS6925444-500V2</SampleSheetName>
  <SampleSheetFolder>D:\Illumina\MiSeq Control Software\SampleSheets</SampleSheetFolder>
  <ManifestFolder>D:\Illumina\MiSeq Control Software\Manifests</ManifestFolder>
  <OutputFolder>D:\Illumina\MiSeqOutput\180823_M03935_0099_000000000-BY7WY</OutputFolder>
  <FocusMethod>AutoFocus</FocusMethod>
  <SurfaceToScan>Both</SurfaceToScan>
  <SaveFocusImages>false</SaveFocusImages>
  <SaveScanImages>true</SaveScanImages>
  <CloudUsername>sebastien.aeby@chuv.ch</CloudUsername>
  <RunManagementType>BaseSpaceCloud</RunManagementType>
  <CloudRunId>115694579</CloudRunId>
  <SendInstrumentHealthToILMN>true</SendInstrumentHealthToILMN>
    '''
    import xml.etree.ElementTree as ET
    root = ET.parse(runinfo_path).getroot()
    run_name = root.findall('ExperimentName')[0].text
    read_length = root.findall('Reads')[0].findall("RunInfoRead")[0].get("NumCycles")
    start_date = datetime.datetime.strptime( root.findall('RunStartDate')[0].text, '%y%m%d').strftime("%Y-%m-%d")

    run2data = {}
    run2data["samples"] = None
    run2data["run_date"] = start_date
    run2data["run_name"] = run_name
    run2data["assay"] = None
    run2data["reads"] = read_length
    run2data["path"] = os.path.dirname(runinfo_path)
    
    return run2data