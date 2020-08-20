

class Jvenn():
    '''
    Generate javascript code for jvenn plots

    - self.add_serie: add a new serie
    - self.add_dictionnatry: add a dictionnary to generate the table displayed under the plot (serie values as keys)
    - self.get_js_code: return complete js code as string

    '''
    
    
    def __init__(self,
                table_header_list,
                href_prefix='',
                href_sufix='',
                table_id="jvenntable",
                venn_div_id="venn_div",
                table_div_id="venn_table"):
        
        '''
        table_header_list: list of values to use as header to the table displayed under the plot.
                           Should match data added with self.add_dictionnatry() 
        href_prefix: hyperlink prefix for each value (optional): eg: chlamdb.ch/ ==> chlamdb.ch/<value>
        href_prefix: hyperlink prefix for each value (optional): eg: /annotation ==> chlamdb.ch/<value>/annotation
        '''
                
        self.dico_string = ''
        
        self.COLOR_LIST = ["rgb(0,102,0)","rgb(90,155,212)","rgb(241,90,96)","rgb(250,220,91)","rgb(255,117,0)","rgb(192,152,83)"]
                
        self.annotation_column = '+ h[this.list[val]]'
        
        self.table_header = f'''<table id="{table_id}" class="table table-striped"><thead class="thead-dark"><tr>%s</tr></thead><tbody>''' % ''.join([f"<th>{i}</th>" for i in table_header_list])
        
        self.jvenn_str = f'''  
var h = new Object();

%s

$(document).ready(function() {{
    $('#{venn_div_id}').jvenn({{
        shortNumber: false,
        colors: %s,
        series: %s,
        displayStat: true,
        displaySwitch: false,
        searchInput:  $("#search-field"),
        searchStatus: $("#search-status"),
        searchMinSize: 1,
        fnClickCallback: function() {{
            var value = '';
            value += '{self.table_header}';
            for (val in this.list) {{
                value += '<tr><td><a href="{href_prefix}' + this.list[val]+ '{href_sufix}' +'">'+ this.list[val]  +'</a></td>' %s +'</tr>';
            }}
            value += '</tbody></table>';
            
            $("#{table_div_id}").html(value);
            
            $('#{table_id}').DataTable( {{
                    dom: 'Bfrtip',
                    "pageLength": 10,
                    "searching": true,
                    "bLengthChange": false,
                    "paging":   true,
                    "info": true
            }} );
    }}
    }});
}});
        '''
        

        self.series = []
    
    
    def add_serie(self,
                  serie_name,
                  value_list):
        '''
        series: [{
           name: "Chlamydia avium 10DC88",
            data: ['group_2561', 'group_2567']
        }, {name: "Chlamydia gallinacea 08-1274/3", data: ['group_2395', 'group_2393']}]
        '''
        serie = {}
        serie["name"] = serie_name
        serie["data"] = value_list
        
        self.series.append(serie)
    
    def add_dictionnatry(self, 
                         dictionnary):
        # //h["group_2395"] = "-</td><td>phospholipase D<br>phospholipase<br>phospholipase D family protein<br>hypothetical protein" ;
        for key, value in dictionnary.items():
            self.dico_string += f'h["{key}"] = "{value}";'

    def get_js_code(self,):
        import re
        
        serie_str = str(self.series)
        
        color_list = self.COLOR_LIST[0:len(self.series)]
        
        if self.dico_string != '':
            code = self.jvenn_str % (self.dico_string,
                                     str(color_list),
                                     serie_str,
                                     self.annotation_column)
        else:
            code = self.jvenn_str % (self.dico_string,
                                     str(color_list),
                                     serie_str,
                                     "")
                                 
        return code