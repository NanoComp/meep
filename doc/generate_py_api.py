#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This tool will extract the docstrings from the meep package and add them to the
hand-written documentation tree.

It is expected that the docstrings will be Markdown compatibile and will not
need to be massaged in any way. One reason for this is to help ensure that the
Look and Feel of the API documentation matches the rest of the documentation.

Before running this script the meep package needs to have been built and the
Python used to run this tool must be able to import meep. Output will be written
to the {project_folder}/doc/docs/Python_User_Interface.md, using the
Python_User_Interface.md.in file as a template.

In order to keep the noise in the module to a minimum, some markdown
snippets/templates used in the generation of the API documentation files can be
found in {project_folder}/doc/_api_snippets, and will be loaded as needed while
processing the docstrings.
"""

import sys
import os
import inspect
import io

import meep

HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
SNIPSDIR = os.path.join(HERE, '_api_snippets')
SRCDOC = os.path.join(HERE, 'docs/Python_User_Interface.md.in')
DESTDOC = os.path.join(HERE, 'docs/Python_User_Interface.md')


# List of names that should not have documentation generated.
#
# Top-level items in the module are automatically excluded if their name is not
# used in a subsition tag in the master template file. Individual class memberes
# can be excluded by using 'Classname.method_name'
EXCLUDES = ['']

#----------------------------------------------------------------------------

class Item(object):
    """
    Common base class for documentable items.
    """

    # Save a copy of the various templates as a class attribute
    templates = dict()

    # This must be overwridden in derived classes
    template_name = None

    def __init__(self, name, obj):
        # load the template if it hasn't been already
        if self.template_name not in Item.templates:
            with open(os.path.join(SNIPSDIR, self.template_name)) as f:
                Item.templates[self.template_name] = f.read()
        self.template = Item.templates[self.template_name]

        # Set other attributes for the item
        self.name = name
        self.obj = obj
        doc = inspect.getdoc(obj)
        if doc:
            doc = inspect.cleandoc(doc)
        self.docstring = doc


    def create_markdown(self, stream):
        raise NotImplementedError


class FunctionItem(Item):
    """
    An introspected item that is a function at module scope.
    """
    template_name = 'function_template.md'

    def __init__(self, name, obj):
        super(FunctionItem, self).__init__(name, obj)
        self.signature = inspect.signature(obj)

    def get_parameters(self, indent):
        sig = inspect.signature(self.obj)
        param_str = str(sig)

        # Wrap and indent the parameters if the line is too long
        if len(param_str) > 50:
            parameters = list(sig.parameters.values())
            params = []
            for idx, param in enumerate(parameters):
                params.append('('+str(param) if idx==0 else ' '*indent+str(param))
            param_str = ',\n'.join(params)
            param_str += ')'
        return param_str

    def create_markdown(self):
        # pull relevant attributes into local variables
        function_name = self.name
        function_name_escaped = function_name.replace('_', '\_')
        docstring = self.docstring if self.docstring else ''
        parameters = self.get_parameters(len(function_name) + 1)

        # Substitute values into the template
        return self.template.format(**locals())


class MethodItem(FunctionItem):
    """
    An introspected item that is a method of a class.
    Mostly the same as a FunctionItem, but can do extra stuff for methods if needed.
    """
    template_name = 'method_template.md'

    def __init__(self, name, obj, klass):
        super(MethodItem, self).__init__(name, obj)
        self.klass = klass
        self.method_name = name
        self.name = '{}.{}'.format(klass.name, name)

    def create_markdown(self):
        # pull relevant attributes into local variables
        class_name = self.klass.name
        method_name = self.method_name
        method_name_escaped = method_name.replace('_', '\_')
        docstring = self.docstring if self.docstring else ''
        parameters = self.get_parameters(4 + len(method_name) + 1)

        # Substitute values into the template
        return self.template.format(**locals())


class ClassItem(Item):
    """
    An introspected item that is a Class.
    """
    template_name = 'class_template.md'

    def __init__(self, name, obj):
        super(ClassItem, self).__init__(name, obj)
        self.add_methods()

    def add_methods(self):
        # Use a match predicate to only look at things in this class, not
        # inherited members
        def _predicate(value):
            if inspect.isfunction(value):
                class_name = value.__qualname__.split('.')[0]
                return class_name == self.name
            else:
                return False

        self.methods = []
        for name, member in inspect.getmembers(self.obj, _predicate):
            if inspect.isfunction(member): # unbound methods are just functions at this point
                self.methods.append(MethodItem(name, member, self))

    def create_markdown(self):
        # pull relevant attributes into local variables
        class_name = self.name
        docstring = self.docstring if self.docstring else ''
        base_classes = [base.__name__ for base in self.obj.__bases__]
        base_classes = ', '.join(base_classes)

        # method_docs = []
        # if self.methods:
        #     # reorder self.methods so __init__ comes first, if it isn't already
        #     methods = self.methods[:]
        #     for idx, meth in enumerate(self.methods):
        #         if meth.name == '__init__':
        #             if idx != 0:
        #                 methods.remove(meth)
        #                 methods.insert(0, meth)
        #             break

        #     for item in methods:
        #         if not check_excluded(item.name) and \
        #            not check_excluded('{}.{}'.format(self.name, item.name)):
        #              doc = item.create_markdown()
        #              method_docs.append(doc)

        # # join the methods into a single string
        # method_docs = '\n'.join(method_docs)

        # Substitute values into the template
        doc = self.template.format(**locals())
        return doc

#----------------------------------------------------------------------------

def check_excluded(name):
    if name in EXCLUDES:
        return True
    if name.startswith('_') and not (name.startswith('__') and name.endswith('__')):
        # It's probably meant to be private
        return True
    return False


def load_module(module):
    """
    Inspect the module and return a list of documentable items in the module.
    """
    items = []
    members = inspect.getmembers(meep)

    for name, member in members:
        if inspect.isclass(member):
            item = ClassItem(name, member)
            items.append(item)
            items += item.methods
        if inspect.isfunction(member):
            items.append(FunctionItem(name, member))

    return items


def generate_docs(items):
    """
    Process the items and create markdown from their docstrings, returned as a
    dictionary.
    """
    docs = dict()
    for item in items:
        if not check_excluded(item.name):
            doc = item.create_markdown()
            docs[item.name] = doc
    return docs


def update_api_document(doc_items):
    """
    Substitute the available class docstrings into the template SRCDOC and write
    to DESTDOC.
    """
    # read
    with open(SRCDOC, 'r') as f:
        srcdoc = f.read()

    # manipulate
    for name, doc in doc_items.items():
        tag = '@@ {} @@'.format(name)
        if tag in srcdoc:
            srcdoc = srcdoc.replace(tag, doc)

    # write results
    with open(DESTDOC, 'w') as f:
        f.write(srcdoc)



def main(args):
    items = load_module(meep)
    doc_items = generate_docs(items)
    update_api_document(doc_items)



if __name__ == '__main__':
    main(sys.argv[1:])