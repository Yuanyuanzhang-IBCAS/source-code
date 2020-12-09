#!/bin/env python
#coding=utf8
'''refRNA main pipeline'''

import os
import sys
import json
import argparse
from time import sleep, asctime
from collections import OrderedDict, defaultdict
from sjmM import execute_cmd, refdir2dict

modulesDir = os.environ['refRNA_dir']+ '/modules'

apps = {
    'QC_v2'             : modulesDir+'/quality_control/qcM.py',
    'Align_v2'          : modulesDir+'/alignment/alignmentM.py',
    'Assembl_v2'        : modulesDir+'/assembl/assemblM.py',
    'AS_v2'             : modulesDir+'/alternative_splicing/asM.py',
    'BamQc_v2'          : modulesDir+'/bam_control/bamM.py',
    'Variation_v2'      : modulesDir+'/call_variation/call_variationM.py',
    'Quantification_v2' : modulesDir+'/quantification/countM.py',
    'Diff_v2'           : modulesDir+'/diff/diffM.py',
    'GO_v2'             : modulesDir+'/go/goM.py',
    'KEGG_v2'           : modulesDir+'/kegg/keggM.py',
    'PPI_v2'            : modulesDir+'/ppi/ppiM.py',
}

def option2conf(opt):
    Basic = OrderedDict()
    if 'projectName' in opt: Basic['projectName'] = opt['projectName']
    if 'taskName' in opt: Basic['taskName'] = opt['taskName']
    if 'taskDir' in opt: Basic['projectDir'] = opt['taskDir']
    if 'queue' in opt: Basic['queue'] = opt['queue']
    if 'taskID' in opt: Basic['taskID'] = opt['taskID']
    if 'species' in opt: Basic['species'] = opt['species']
    Basic['modules']=opt['mod'] if 'mod' in opt else apps.keys() 
    #Basic['modules']=apps.keys()
    if 'seqType' in opt: Basic['seqType'] = opt['seqType']
    if 'ssLib' in opt: Basic['ssLib'] = opt['ssLib']
    if 'dataDir' in opt: Basic['rawDir'] = opt['dataDir']
    if 'samples' in opt: 
        ds = {}
        for ss in opt['samples']:
            fn_lst = []
            fn = ss.keys()[0] 
            for s in ss.values()[0]:
                if s != None:
                    fn_lst.append(Basic['rawDir']+'/'+s)
            ds[fn] = fn_lst
        Basic['fq'] = ds
                    
        #Basic['fq'] = { s.keys()[0]:[ Basic['rawDir']+'/'+f for f in s.values()[0]] for s in opt['samples']}
    if 'groups' in opt: Basic['groups'] = opt['groups']
    if 'compares' in opt: Basic['compares'] = opt['compares']
    if 'venn' in opt: Basic['venn'] = opt['venn']
    
    if 'species' in opt:
        Basic['Reference']=OrderedDict
        Basic['Reference']=refdir2dict(refdir = opt['species'])
        if 'kegg' in opt: Basic['Reference']['kegg'] = opt['kegg']
        if 'ppi' in opt: Basic['Reference']['ppi'] = opt['ppi']
    
    ToolsDict=OrderedDict()
    if 'alignTool' in opt: ToolsDict['alignTool'] = opt['alignTool']
    if 'phred' in opt: ToolsDict['phred'] = opt['phred']
    if 'assemlTool' in opt: ToolsDict['assemlTool'] = opt['assemlTool']
    if 'diffTool' in opt: ToolsDict['diffTool'] = opt['diffTool']
    if 'foldchange' in opt: ToolsDict['foldchange'] = opt['foldchange']
    if 'pType' in opt: ToolsDict['pType'] = opt['pType']
    if 'pValue' in opt: ToolsDict['pValue'] = opt['pValue']
    if ToolsDict: Basic['Tools'] = ToolsDict
    
    Basic['submit_time'] = asctime() 

    paras_dict = OrderedDict()
    paras_dict['Basic'] = Basic
    
    with open(Basic['projectDir']+'/config.json', 'w') as fh:
        json.dump(paras_dict, fh, indent=4)
    return Basic['projectDir']

def module_order(fileName):
    #order={ i['nodeName']: [j['name'] for j in i['linkData']] for i in opt['order']}
    order = json.load(open(fileName, 'r'))['order']
    new_order = defaultdict(list)
    for i in order:
        for j in order[i]:
            new_order[j].append(i)
    return new_order

def reportData_ini(paras, logDir):
    # for report
    render = OrderedDict({
        'nodeDataArray':'''[
        {key: "基本质控", loc:"0 0", color:finishedColor},
        {key: "序列比对", loc:"-1 100", color:finishedColor},
        {key: "高级质控", loc:"150 100", color:finishedColor},
        {key: "SNP和InDel检测", loc:"123 200", color:finishedColor},
        {key: "转录本组装", loc:"-7 200", color:finishedColor},
        {key: "可变剪接分析", loc:"-150 200", color:finishedColor},
        {key: "差异表达分析", loc:"-15 300", color:finishedColor},
        {key: "蛋白质互作分析", loc:"-165 400", color:finishedColor},
        {key: "GO富集分析", loc:"-13 400", color:finishedColor},
        {key: "KEGG富集分析", loc:"135 400", color:finishedColor}
        ]''',

        'linkDataArray':'''[
        {from: "基本质控", to:"序列比对", "curviness":0},
        {from: "序列比对", to:"高级质控", "curviness":0},
        {from: "高级质控", to:"SNP和InDel检测", "curviness":0},
        {from: "序列比对", to:"转录本组装", "curviness":0},
        {from: "转录本组装", to: "可变剪接分析", "curviness":0},
        {from: "转录本组装", to: "差异表达分析", "curviness":0},
        {from: "差异表达分析", to: "蛋白质互作分析", "curviness":-20},
        {from: "差异表达分析", to: "GO富集分析", "curviness":0},
        {from: "差异表达分析", to: "KEGG富集分析", "curviness":20}
        ]'''
    })    

    render['flag_paras']=True
    render.update(paras)

    with open(logDir+'/.reportData.json','w') as fh:
        json.dump(render, fh, indent=4)
 
    
if __name__== "__main__":
    parser = argparse.ArgumentParser(description="refRNA main pipeline ")
    parser.add_argument('--param_val', type=json.loads, help='json', required=True)
    args = vars(parser.parse_args())

    taskDir=option2conf(args['param_val'])

    order_json = os.environ['refRNA_dir']+ '/order.json'
    order=module_order(order_json)

    config=taskDir+'/config.json'

    paras = json.load(open(config,'r'), object_pairs_hook=OrderedDict)
    #taskDir = paras['Basic']['taskDir']
    logDir = taskDir + '/Log'
    taskID = paras['Basic']['taskID']
    
    assert not os.system('mkdir -p %s'%(logDir))
    
    # 参数数据
    reportData_ini(paras['Basic'], logDir)

    # generate sjm
    sjmfiles=[]
    for m in paras['Basic']['modules']:
        sjmFile=logDir + '/.%s_sjm.job'%(m)
        cmd='python {app} --anamod {mod} --config {config} --sjmFile {sjmFile}'.format(app=apps[m], mod=m, config=config, sjmFile=sjmFile)
        execute_cmd(cmd)
        sjmfiles.append(sjmFile)
    
    libDir = os.environ['refRNA_dir']+ '/lib'
    outjob = logDir + '/final_sjm.job'
    cmd="perl {app} --jobs {jobfiles} --orders '{order}' --out {outjob}".format(app=libDir+'/cat_sjms.pl',
        jobfiles=','.join(sjmfiles), order=json.dumps(order), outjob=outjob)  
    execute_cmd(cmd)

    '''
    if 'queue' in  paras['Basic']:
        cmd = "sed -i 's/-q all.q/{new_queue}/' {sjm_file}".format(new_queue=paras['Basic']['queue'], sjm_file=outjob)
        execute_cmd(cmd)
        sleep(1)
    '''
    #cmd = "sed -i 's/-q all.q/{new_queue}/' {sjm_file}".format(new_queue='-q login04.q', sjm_file=outjob)
    #execute_cmd(cmd)
    #sleep(1)

    
    assert not os.system('sjm --save_interval 5 --check_interval 5 -x %s  %s'%(taskID, outjob))

