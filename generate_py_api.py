"""
This tool will extract the docstrings from the meep package and add them to the
hand-written documentation tree.

It is expected that the docstrings will be Markdown compatibile and will not
need to be massaged in any way. One reason for this is to help ensure that the
Look and Feel of the API documentation matches the rest of the documentation.

Before running this script the meep package needs to have been built and the
Python used to run this tool must be able to import meep. Output will be written
to the {project_folder}/doc/docs/Python_User_Interface.md, using the
Python_User_Interface.md.in file asa template.

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

SNIPSDIR = 'doc/_api_snippets'
SRCDOC = 'doc/docs/Python_User_Interface.md.in'
DESTDOC = 'doc/docs/Python_User_Interface.md'

# List of names that should not have documentation generated.
#
# Top-level items in the module are automatically excluded if their name is not
# used in a subsition tag in the master template file. Individual class memberes
# can be excluded by using 'Classname.Methodname'
EXCLUDES = ['']

#----------------------------------------------------------------------------

class Item(object):
    def __init__(self, name, obj):
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
    def __init__(self, name, obj):
        super().__init__(name, obj)
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

    def create_markdown(self, stream):
        pass


class MethodItem(FunctionItem):
    """
    An introspected item that is a method of a class.
    Mostly the same as a FunctionItem, but can do extra stuff for methods if needed.
    """
    def __init__(self, name, obj, klass):
        super().__init__(name, obj)
        self.klass = klass
        with open(os.path.join(SNIPSDIR, 'method_template.md')) as f:
            self.template = f.read()

    def create_markdown(self):
        # pull relevant attributes into local variables
        class_name = self.klass.name
        method_name = self.name
        method_name_escaped = method_name.replace('_', '\_')
        docstring = self.docstring if self.docstring else ''
        parameters = self.get_parameters(4 + len(method_name) + 1)

        # Substitute values into the template
        return self.template.format(**locals())



class ClassItem(Item):
    """
    An introspected item that is a Class.
    """
    def __init__(self, name, obj):
        super().__init__(name, obj)
        self.add_methods()
        with open(os.path.join(SNIPSDIR, 'class_template.md')) as f:
            self.template = f.read()

    def add_methods(self):
        self.methods = []
        for name, member in inspect.getmembers(self.obj):
            if inspect.isfunction(member):
                self.methods.append(MethodItem(name, member, self))

    def create_markdown(self):
        # pull relevant attributes into local variables
        if self.name == 'Simulation':
            print('break here')

        class_name = self.name
        docstring = self.docstring if self.docstring else ''
        base_classes = [base.__name__ for base in self.obj.__bases__]
        base_classes = ', '.join(base_classes)

        method_docs = []
        if self.methods:
            for item in self.methods:
                if not check_excluded(item.name) and \
                   not check_excluded(f'{self.name}.{item.name}'):
                     doc = item.create_markdown()
                     method_docs.append(doc)
        method_docs = '\n'.join(method_docs)

        # TODO: reorder self.methods so __init__ comes first

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
            items.append(ClassItem(name, member))
            if name == 'Simulation':
                print('break here')
        if inspect.isfunction(member):
            items.append(FunctionItem(name, member))

    return items


def generate_class_docs(items):
    """
    Process the items (just classes for now) and create markdown from their
    docstrings, returned as a dictionary.
    """
    class_docs = dict()
    for item in items:
        if isinstance(item, ClassItem) and not check_excluded(item.name):
            doc = item.create_markdown()
            class_docs[item.name] = doc
    return class_docs


def update_api_document(class_items):
    """
    Substitute the available class docstrings into the template SRCDOC and write
    to DESTDOC.
    """
    # read
    with open(SRCDOC, 'r') as f:
        srcdoc = f.read()

    # manipulate
    for name, doc in class_items.items():
        tag = f'@@ {name} @@'
        if tag in srcdoc:
            srcdoc = srcdoc.replace(tag, doc)

    # write results
    with open(DESTDOC, 'w') as f:
        f.write(srcdoc)


def main(args):
    items = load_module(meep)

    class_items = generate_class_docs(items)
    update_api_document(class_items)



if __name__ == '__main__':
    main(sys.argv[1:])