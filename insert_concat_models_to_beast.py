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

def insert_patterns_for_locus(xml, models, insert_position, model_data):
    sites = model_data.find('partition')[0]
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

def append_to_element_name(element, key, text):
    return '.'.join([text, element.get(key)])


def strip(element, text, tail = None, comment = False):
    if not tail:
        tail = text
    if element.text and not comment:
        l = len(element.text.split("\t"))
        element.text = "\n" + "\t" * (l - text)
    if element.tail:
        l = len(element.tail.split("\t"))
        element.tail = "\n" + "\t" * (l - tail)
    return element

def iterate_over_model_children(element, parent, model, locus, site = None, dedent = 2):
    #pdb.set_trace()
    try:
        if element.tag.func_name == 'Comment':
            element = strip(element, dedent + 1, None, True)
    except:
        if element.get('id') == model.upper():
            name = append_to_element_name(element, 'id', locus)
            element.set('id', name)
            element = strip(element, dedent)
        elif element.get('idref'):
            #pdb.set_trace()
            name = append_to_element_name(element, 'idref', locus)
            element.set('idref', name)
            element = strip(element, dedent + 1)
        elif element.tag == 'parameter':
            #
            name = append_to_element_name(element, 'id', locus)
            element.set('id', name)
            element = strip(element, dedent + 1)
        else:
            element = strip(element, dedent + 1)
    for child in element:
        iterate_over_model_children(child, element, model, locus)
    return element

def insert_comments(xml, insert_position, text):
    comment = Comment(text="{0} automatically inserted by beast_tools".format(text))
    comment.tail = '\n\n\t'
    xml.insert(insert_position, comment)
    insert_position += 1
    return xml, insert_position

def get_subs_model_name(model):
    m = model.split('-')[1]
    return m.rstrip('IG').lower()

def get_site_model_name(model):
    return model.split('-')[1].lower()

def insert_operators_for_locus(xml, models, snippets):
    operators = xml.find('operators')
    c = 0
    for name in models:
        m_name = get_subs_model_name(name)
        s_name = get_site_model_name(name)
        model_operators = snippets.find(m_name).find('operators')
        for locus in models[name]:
            frequencies = model_operators.find('frequencies')
            if frequencies:
                frequencies_copy = copy.deepcopy(frequencies)
                frequencies_copy = iterate_over_model_children(frequencies_copy, None, m_name, locus, 0)
                pdb.set_trace()
            transtrans = model_operators.find('transtrans')
            if transtrans:
                transtrans_copy = copy.deepcopy(transtrans)
                transtrans_copy = iterate_over_model_children(transtrans_copy, None, m_name, locus, 0)
            operators.insert(-1, frequencies_copy)
            operators.insert(-1, transtrans_copy)
    return xml

def insert_models_for_locus(xml, models, insert_position, snippets):
    xml, insert_position = insert_comments(xml, insert_position, "Begin models")
    for m in models:
        m_name = get_subs_model_name(m)
        model = snippets.find(m_name)
        structure = model.find('model')[0]
        for locus in models[m]:
            # create a new object copy to edit on every iteration
            structure_copy = copy.deepcopy(structure)
            structure_copy = iterate_over_model_children(structure_copy, None, m_name, locus)
            xml.insert(insert_position, structure_copy)
            insert_position += 1
            #pdb.set_trace()
    xml, insert_position = insert_comments(xml, insert_position, "End models")
    return xml, insert_position

def insert_site_models_for_locus(xml, models, insert_position, snippets):
    xml, insert_position = insert_comments(xml, insert_position, "Begin sites")
    for name in models:
        m_name = get_subs_model_name(name)
        s_name = get_site_model_name(name)
        site_model = snippets.find(m_name).find('sites').find(s_name)[0]
        for locus in models[name]:
            site_model_copy = copy.deepcopy(site_model)
            site_model_copy = iterate_over_model_children(site_model_copy, None, s_name, locus, m_name, 3)
            xml.insert(insert_position, site_model_copy)
            insert_position += 1
    xml, insert_position = insert_comments(xml, insert_position, "End sites")
    return xml, insert_position

def iterate_over_likelihood(element, m_name, s_name, locus):
    element.set('id', "{0}.likelihood".format(locus))
    for child in element:
        if child.tag == 'patterns':
            child.set('idref', locus)
        #elif child.tag == 'treeModel':
        #    child.set('idref', "{0}.{1}".format(locus, m_name.upper()))
        elif child.tag == 'siteModel':
            child.set('idref', "{0}.{1}".format(locus, s_name.upper()))
        else:
            pass
    return element

def insert_tree_likelihoods_for_locus(xml, models, insert_position, snippets):
    xml, insert_position = insert_comments(xml, insert_position, "Begin likelihoods")
    for name in models:
        m_name = get_subs_model_name(name)
        s_name = get_site_model_name(name)
        likelihood = snippets.find('likelihood')[0]
        for locus in models[name]:
            likelihood_copy = copy.deepcopy(likelihood)
            likelihood_copy = iterate_over_likelihood(likelihood_copy, m_name, s_name, locus)
            xml.insert(insert_position, likelihood_copy)
            insert_position += 1
    xml, insert_position = insert_comments(xml, insert_position, "End likelihoods")
    return xml, insert_position

def get_items_to_recursively_delete(element, value, to_delete, parent = None):
    if element.tag == 'parameter' and element.get('idref') in value:
        to_delete.append(parent)
    for child in element:
        get_items_to_recursively_delete(child, value, to_delete, element)
    return to_delete

def delete_children_from_node(xml, section, value):
    to_delete = get_items_to_recursively_delete(section, value, [])
    [section.remove(d) for d in to_delete]
    return xml

def get_file_log_section(xml):
    log_sections = xml.find('mcmc').findall('log')
    for section in log_sections:
        if section.get('id') == 'fileLog':
            return section[0]
    
def main():
    super_concat = False
    options, args = interface()
    xml = ElementTree().parse(options.input, parser=XMLParser(target=MyTreeBuilder()))
    # delete the older subs. models from the xml file
    for node in ['HKYModel', 'gtrModel','siteModel', 'patterns', 'treeLikelihood']:
        xml = delete_node(xml, node, 1) 
    # delete the kappa and frequency parameters in 'operators'
    delete = ['kappa', 'frequencies', 'alpha', 'pInv', 'ac', 'ag', 'at', 'cg', 'gt']
    xml = delete_children_from_node(xml, xml.find('operators'), delete)
    xml = delete_children_from_node(xml, xml.find('mcmc').find('posterior').find('prior'), delete)
    # there are 2 log tags, disambiguated w/ id params.  delete elements from 
    # the one we want (fileLog)
    xml = delete_children_from_node(xml, get_file_log_section(xml), delete)
    
    # jettison some comments
    xml = comment_remover(xml, ['The general time reversible', 'HKY substitution model','site model', 'The unique patterns from 1 to end', 'npatterns=']) 
    
    # load our substitution model information
    substitution_models = cPickle.load(open(options.subs))
    snippets = ElementTree().parse(options.params, parser=XMLParser(target=MyTreeBuilder()))
    
    # insert patterns on a per locus basis
    insert_position = get_position(xml, 'alignment')
    xml = insert_patterns_for_locus(xml, substitution_models, insert_position, snippets)
    
    # insert substitution models on a per locus basis
    insert_position = get_position(xml, 'strictClockBranchRates')
    xml, insert_position = insert_models_for_locus(xml, substitution_models, insert_position, snippets)
    
    # insert site models on a per locus basis
    xml, insert_position = insert_site_models_for_locus(xml, substitution_models, insert_position, snippets)
    
    # insert tree likelihoods on a per locus basis
    xml, insert_position = insert_tree_likelihoods_for_locus(xml, substitution_models, insert_position, snippets)
    
    # insert operators
    xml = insert_operators_for_locus(xml, substitution_models, snippets)
    write(xml, options.output)
    pdb.set_trace()
    
    """model_names, site_names = get_xml_model_names(set(sub_models_from_modeltest.values()))
    
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
    """
    


if __name__ == '__main__':
    main()

