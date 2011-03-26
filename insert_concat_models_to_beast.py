#!/usr/bin/env python
# encoding: utf-8
"""
insert_models_to_beast.py

Created by Brant Faircloth on 18 February 2011 14:15 PST (-0800).
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import copy
import cPickle
import optparse
from xml.etree.cElementTree import ElementTree
from xml.etree.ElementTree import TreeBuilder
from xml.etree.ElementTree import Comment
from xml.etree.cElementTree import XMLParser
from xml.etree.cElementTree import tostring
from StringIO import StringIO


import pdb

def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)
   
    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to the input file.', 
metavar='FILE')

    p.add_option('--output', dest = 'output', action='store', 
type='string', default = None, help='The path to the input file.', 
metavar='FILE')

    p.add_option('--new-params', dest = 'params', action='store', 
type='string', default = None, help='The path to the params xml file.', 
metavar='FILE')

    p.add_option('--sub-models', dest = 'subs', action='store', 
type='string', default = None, help='The path to the substitution models file.', 
metavar='FILE')

    (options, args) = p.parse_args()
    for k,v in options.__dict__.items():
        if not v:
            print "{0} is missing a value".format(k)
            p.print_help()
            sys.exit(2)
    return options, args

class MyTreeBuilder(TreeBuilder):
   def comment(self, data):
       self.start(Comment, {})
       self.data(data)
       self.end(Comment)

def get_level(level, node):
    level = "." + '/' * level
    level = '{0}{1}'.format(level, node)
    return level

def delete_node(xml, node, level = 1):
    level = get_level(level, node)
    for node in xml.findall(level):
        xml.remove(node)
    return xml

def delete_children_from_node(xml, node, child_parameter, level=1):
    level = get_level(level, node)
    search = xml.find(level)
    for node in search:
        for child in iter(node):
            if child.tag == 'parameter' and child.get('idref') == child_parameter:
                search.remove(node)
    return xml

def delete_children_from_log_node(xml, node, child_parameter, level = 2):
    level = get_level(level, node)
    search = xml.findall(level)
    #pdb.set_trace()
    for node in search:
        if node.get('id') == 'fileLog':
            for child in iter(node):
                if child.tag == 'parameter' and child.get('idref') == child_parameter:
                    node.remove(child)
    return xml

def model_name_formatter(m):
    if len(m) == 3 or m == 'JC69':
        return 'siteModel_{0}'.format(m)
    # why the F does this one model have a 4-letter abbr.
    elif m == 'JC69I' or m == 'JC69G':
        return 'siteModel_{0}_{1}'.format(m[0:4], m[4])
    elif len(m) == 4:
        return 'siteModel_{0}_{1}'.format(m[0:3], m[3])
    elif m == 'JC69IG':
        return 'siteModel_{0}_{1}_{2}'.format(m[0:4], m[4], m[5])
    elif len(m) == 5:
        return 'siteModel_{0}_{1}_{2}'.format(m[0:3], m[3], m[4])

def get_xml_model_names(models):
    xml_model_names = set([])
    xml_site_model_names = set([])
    for m in models:
        if m == 'JC69':
            xml_model_names.add(m)
        else:
            xml_model_names.add(m[0:3])
        xml_site_model_names.add(model_name_formatter(m))
    return xml_model_names, xml_site_model_names

def get_generic_section_to_add(model_names, model_data, section):
    """get the specific model data we are going to add to the beat file"""
    models_to_add = []
    models = model_data.find(section)
    for m in iter(models):
        if m.get('id') in model_names:
            m.text = m.text.replace('\t\t\t','\t')
            m.tail = m.tail.replace('\t\t','\t')
            models_to_add.append(m)
    return models_to_add

def get_generic_section_children_to_add(model_names, model_data, section):
    """get the model-specific operators to insert"""
    #if section == 'operators':
    #    pdb.set_trace()
    operators_to_add = []
    operators = model_data.find(section)
    for node in iter(operators):
        for child in iter(node):
            if child.tag == 'parameter' and child.get('idref').split('.')[0] in model_names:
                operators_to_add.append(node)
    return operators_to_add

def get_log_entries_to_add(model_names, model_data):
    """docstring for log_entries_to_add"""
    log_entries_to_add = []
    log_entries = model_data.find('log')
    for node in iter(log_entries):
        if node.tag == 'parameter' and node.get('idref').split('.')[0] in model_names:
            log_entries_to_add.append(node)
    return log_entries_to_add

def comment_remover(xml, comments):
    for comment in comments:
        for element in iter(xml):
            try:
                if element.tag.func_name == 'Comment':
                    if comment in element.text:
                        xml.remove(element)
            except AttributeError:
                pass
    return xml

def write(xml, outf):
    """use pretty print because other options suck"""
    f = open(outf, 'w')
    f.write(tostring(xml))
    f.close()

def get_position(xml, tag):
    for k, child in enumerate(xml.getchildren()):
        if child.tag == tag:
            insert_position = k + 1
    return insert_position

def insert_models_and_sites(xml, insert_position, models, sites):
    xml.insert(insert_position, Comment(text="Begin model and site values automatically"+
                " inserted by insertobeast.py"))
    insert_position += 1
    for model in models:
        xml.insert(insert_position, model)
        insert_position += 1
    xml.insert(insert_position, Comment(text="Site values automatically"
                +" inserted by insertobeast.py"))
    insert_position += 1
    for site in sites:
        xml.insert(insert_position, site)
        insert_position += 1
    xml.insert(insert_position, Comment(text="End model and site values automatically"
                +" inserted by insertobeast.py"))
    insert_position += 1
    return insert_position

def insert_tree_likelihoods(xml, insert_position, models):
    xml.insert(insert_position, Comment(text="Begin tree likelihoods automatically"+
                " inserted by insertobeast.py"))
    insert_position += 1
    #for model

def update_tree_likelihoods(xml, sub_models, level = 1):
    """docstring for update_tree_likelihoods"""
    # convert dict
    locus_dict = {}
    for k,v in sub_models.iteritems():
        locus_dict[k] = model_name_formatter(v)
    level = get_level(level, 'treeLikelihood')
    search = xml.findall(level)
    for node in search:
        for child in iter(node):
            if child.tag == 'siteModel':
                #pdb.set_trace()
                locus = node.get('id').split('.')[0]
                child.set('idref', locus_dict[locus])
    return xml

def insert_to_generic_sections(xml, values, section, id_name, level = 2):
    begin_comm = Comment(text="Begin {0} values automatically inserted by insertobeast.py".format(section))
    end_comm = Comment(text="End {0} values automatically inserted by insertobeast.py".format(section))
    values.insert(0, begin_comm)
    values.append(end_comm)
    level = get_level(level, section)
    search = xml.findall(level)
    for result in search:
        if result.get('id') == id_name:
            for k, v in enumerate(values):
                result.insert(k, v)
    return xml

def insert_sites_from_models(xml, models, insert_position, model_data):
    sites = model_data.find('patternSection')[0]
    test = []
    for model in models:
        for locus in models[model]:
            site_copy = copy.deepcopy(sites)
            site_copy.set('id', locus)
            site_copy.set('from', str(models[model][locus][0]))
            site_copy.set('to', str(models[model][locus][1]))
            xml.insert(insert_position, site_copy)
            insert_position += 1
    return xml

def main():
    super_concat = False
    options, args = interface()
    xml = ElementTree().parse(options.input, parser=XMLParser(target=MyTreeBuilder()))
    # delete the older subs. models from the xml file
    for node in ['HKYModel','siteModel', 'patterns']:
        xml = delete_node(xml, node, 1) 
    if super_concat:
        xml = delete_node(xml, 'treeLikehood', 1)
    # delete the kappa and frequency parameters in 'operators'
    for parameter in ['kappa', 'frequencies']:
        xml = delete_children_from_node(xml, 'operators', parameter)
        xml = delete_children_from_node(xml, 'prior', parameter, 2)
        xml = delete_children_from_log_node(xml, 'log', parameter)
    # jettison some comments
    xml = comment_remover(xml, ['HKY substitution model','site model', 'The unique patterns from 1 to end', 'npatterns=']) 
    # get our subs model information
    #sub_models_from_modeltest = {line.strip().split('\t')[0]:line.strip().split('\t')[1].split('-')[1]
    #                                for line in open(options.subs, 'rU')}
    sub_models_from_modeltest = cPickle.load(open(options.subs))
    model_data = ElementTree().parse(options.params, parser=XMLParser(target=MyTreeBuilder()))
    
    # generate sites xml
    insert_position = get_position(xml, 'alignment')
    xml = insert_sites_from_models(xml, sub_models_from_modeltest, insert_position, model_data)
    write(xml, options.output)
    pdb.set_trace()
    model_names, site_names = get_xml_model_names(set(sub_models_from_modeltest.values()))
    
    # get the xml data that we need to add for the models and their parameters
    models_to_add = get_generic_section_to_add(model_names, model_data, 'models')
    sites_to_add = get_generic_section_to_add(site_names, model_data, 'sites')
    operators_to_add = get_generic_section_children_to_add(model_names, model_data, 'operators')
    log_entries_to_add = get_log_entries_to_add(model_names.union(site_names), model_data)
    priors_to_add = get_generic_section_children_to_add(model_names, model_data, 'priors')
    
    # get the last position of the strictClockBranchRates
    insert_position = get_position(xml, 'strictClockBranchRates')
    # insert the models and sites we need
    insert_position = insert_models_and_sites(xml, insert_position, models_to_add, sites_to_add)
    # modify the tree likelihood statements
    if not super_concat:
        xml = update_tree_likelihoods(xml, sub_models_from_modeltest)
    else:
        insert_position = get_position(xml, 'siteModel')
        xml = insert_tree_likelihoods(xml, sub_models_from_modeltest)
    # insert the operators we need
    xml = insert_to_generic_sections(xml, operators_to_add, 'operators', 'operators')
    # insert the priors we need
    xml = insert_to_generic_sections(xml, priors_to_add, 'prior', 'prior')
    # alter the log node to collect data
    xml = insert_to_generic_sections(xml, log_entries_to_add, 'log', 'fileLog')
    # write to the output file
    


if __name__ == '__main__':
    main()

