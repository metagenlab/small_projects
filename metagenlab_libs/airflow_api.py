
import requests
from requests.auth import HTTPBasicAuth

class AirflowApi():
    
    
    def __init__(self, username, password):
        
        self.API_ENDPOINT = "http://0.0.0.0:8081/api/experimental/"
        
        self.headers = {'Content-Type': 'application/json',
                        'Cache-Control': 'no-cache'}
        
        self.username = username 
        
        self.password = password 
        
        
        
    def run_dag(self, 
                dag_id, 
                dag_json=False):
        
        t = f'dags/{dag_id}/dag_runs'
        
        if dag_json:
            r = requests.post(url = self.API_ENDPOINT + t, json = dag_json, auth=HTTPBasicAuth(self.username, self.password), headers=self.headers) 
        else:
            r = requests.post(url = self.API_ENDPOINT + t, auth=HTTPBasicAuth(self.username, self.password), headers=self.headers) 
        
        # {'execution_date': '2020-11-19T16:10:56+00:00', 'message': 'Created <DagRun epidemiology_with_jbrowse @ 2020-11-19 16:10:56+00:00: manual__2020-11-19T16:10:56+00:00, externally triggered: True>', 'run_id': 'manual__2020-11-19T16:10:56+00:00'}    
        res = eval(r.text)
        
        return res
        
        
    def check_status(self, 
                     dag_id,
                     run_id):
        
        print("run id",run_id)
        run_date = run_id.split("__")[1]
        
        t = f'dags/{dag_id}/dag_runs/{run_date}'
        r = requests.get(url = self.API_ENDPOINT + t, auth=HTTPBasicAuth(self.username, self.password))
        res = eval(r.text)

        return res["state"]
    
    
    def retrieve_all_run_data(self, 
                              dag_id):
        
        t = f'dags/{dag_id}/dag_runs'
        r = requests.get(url = self.API_ENDPOINT + t, auth=HTTPBasicAuth(self.username, self.password))
        res = eval(r.text)

        return res
    
    def retrieve_run_data(self,
                          dag_id,
                          run_id):
               
        runs_data = self.retrieve_all_run_data(dag_id)
        
        target_run = [i for i in runs_data if i["run_id"]==run_id]
        
        if len(target_run) > 1:
            raise IOError("Unexpected number of hits")

        return target_run[0]
        