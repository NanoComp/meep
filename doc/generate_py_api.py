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

The following tag patterns are used to indicate where in the input template that
docstrings will be inserted:

  - @@ top_level_function_name @@
  - @@ ClassName @@
  - @@ ClassName.method_name @@
  - @@ ClassName[all-methods] @@
  - @@ ClassName[methods-with-docstrings] @@

If a docstring containes text that should be treated as an alternate function
signature, then those lines can be marked with tags like `##sig` (to move the
line to the header) or `##sig-keep` to copy the line, but not remove it from the
docstring, keeping it in place.

    NOTE: Currently the signature tags functionality assumes that the signature
          does not span more than one line.

In order to get the actual text for default parameters, (rather than the value
that the inspect module gives us) Python's ast module is used to find the actual
text from the source code for those defaults. Then means that we can get
constants like `mp.ALL_COMPONENTS` instead of `20`. This works well for most
cases, but the AST for the following types of default parameters are complex enough
that it's probably not worth the effort to use the AST for them. In these cases
they fallback to using the `inspect` module's evaluated form of the default.
With a little tweaking of the code the downsides of this approach can be
mitigated.

* Lambdas: The easiest way to workaround this is to assign the lambda to a
  global and then use that global as the default value for the parameter.

* Function calls: If the object returned by the call has a good `__repr__`
  representation, then that will be used for the default's text in the
  documentaion. If not, then an approach like the above is recommended.

* Math expressions: An example of this case is `angle=2*math.pi` which will
  result in the documentation using `angle=6.283185307179586`. If it's desired
  to have something more meaningful in the API documentation then assigning the
  value to a symbolic name and using that name for the default value works well
  in this case too.
"""

import sys
import os
import inspect
import io
import ast
import itertools
import textwrap

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


# ast.Constant was added in Python 3.6, and ast.NameConstant is deprecated
# starting in Python 3.8, so use Constant if it exists, otherwise fall back to
# NameConstant.
try:
    Constant = ast.Constant
except AttributeError:
    Constant = ast.NameConstant

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
        self.sig = inspect.signature(self.obj)
        try:
            # Try using the AST to get the actual text for default values
            parameters = self.get_parameters_from_ast()
            param_str = '({})'.format(', '.join(parameters))

            # Wrap and indent the parameters if the line is too long
            if len(param_str) > 50:
                params = []
                for idx, param in enumerate(parameters):
                    params.append('('+str(param) if idx==0 else ' '*indent+param)
                param_str = ',\n'.join(params)
                param_str += ')'
            return param_str

        except ValueError as ex:
            print("Warning: falling back to old parameter extraction method for {}\n{}".format(self.name, ex))
        except OSError:
            # This happens when there is no source for some function/class (like NamedTuples)
            pass

        # Fall back to using just the inspect module's info
        param_str = str(self.sig)

        # Wrap and indent the parameters if the line is too long
        if len(param_str) > 50:
            parameters = list(self.sig.parameters.values())
            params = []
            for idx, param in enumerate(parameters):
                params.append('('+str(param) if idx==0 else ' '*indent+str(param))
            param_str = ',\n'.join(params)
            param_str += ')'
        return param_str


    def get_parameters_from_ast(self):
        # Use Python's ast module to parse the function's code and pull out parameter names and
        # the text for default values.
        def parse_name(node):
            assert isinstance(node, ast.arg)
            if node.annotation != None:
                raise ValueError("Annotations are not currently supported")
            return node.arg

        source = inspect.getsource(self.obj)
        source = textwrap.dedent(source)
        module = ast.parse(source)
        func = module.body[0]

        if not isinstance(func, ast.FunctionDef):
            raise ValueError('Unsupported node type: {}'.format(type(func)))

        # Get the param name and defaults together into a list of tuples
        args = reversed(func.args.args)
        defaults = reversed(func.args.defaults)
        iter = itertools.zip_longest(args, defaults, fillvalue=None)
        parameters = [(parse_name(name), default) for name, default in reversed(list(iter))]

        # Convert each of the parameters (name, ast.node) to valid Python parameter
        # code, like what would have been in the actual source code.
        transformed = []
        for name, node in parameters:
            if node is None:
                transformed.append(name)
            else:
                default = self.transform_node(name, node)
                transformed.append("{}={}".format(name, default))

        # Check for vararg (like *foo) and kwarg (like **bar) parameters
        if func.args.vararg is not None:
            transformed.append('*{}'.format(func.args.vararg.arg))
        if func.args.kwarg is not None:
            transformed.append('**{}'.format(func.args.kwarg.arg))

        return transformed

    def transform_node(self, name, node):
        if isinstance(node, (ast.Str, ast.Bytes)):
            default = repr(node.s)
        elif isinstance(node, ast.Num):
            default = repr(node.n)
        elif isinstance(node, Constant):
            default = node.value
        elif isinstance(node, ast.Name):
            default = node.id
        elif isinstance(node, ast.Attribute):
            a = []
            n = node
            while isinstance(n, ast.Attribute):
                a.append(n.attr)
                n = n.value
            assert isinstance(n, ast.Name), "I'm confused..."
            a.append(n.id)
            default = ".".join(reversed(a))
        # elif isinstance(node, ast.Lambda):
        #     default = 'FIXME_Lambda'
        # elif isinstance(node, ast.Call):
        #     default = 'FIXME_Call'
        else:
            param = self.sig.parameters[name]
            default = param.default
            # print('Using inspect module for {}={} in {}'.format(name, default, self.name))
        return default


    def check_other_signatures(self, docstring):
        """ Search for alternate function signatures in the docstring """
        SIG = '##sig'
        SIG_KEEP = '##sig-keep'
        other_signatures = ''

        if SIG in docstring or SIG_KEEP in docstring:
            other_sigs = []
            docstring_lines = []
            for line in docstring.split('\n'):
                if line.endswith(SIG):
                    line = line.replace(SIG, '')
                    line = line.strip()
                    line = line.replace('`', '')
                    other_sigs.append(line)

                elif line.endswith(SIG_KEEP):
                    line = line.replace(SIG_KEEP, '')
                    line = line.strip()
                    docstring_lines.append(line)
                    line = line.replace('`', '')
                    other_sigs.append(line)

                else:
                    docstring_lines.append(line)

            docstring = '\n'.join(docstring_lines)
            other_signatures = ''.join(['\ndef {}:'.format(item) for item in other_sigs])
        return docstring, other_signatures

    def create_markdown(self):
        # pull relevant attributes into local variables
        function_name = self.name
        function_name_escaped = function_name.replace('_', '\\_')
        docstring = self.docstring if self.docstring else ''
        parameters = self.get_parameters(len(function_name) + 1)
        docstring, other_signatures = self.check_other_signatures(docstring)

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
        method_name_escaped = method_name.replace('_', '\\_')
        docstring = self.docstring if self.docstring else ''
        parameters = self.get_parameters(4 + len(method_name) + 1)
        docstring, other_signatures = self.check_other_signatures(docstring)

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

        # Substitute values into the template
        docs = dict()
        class_doc = self.template.format(**locals())
        docs[class_name] = class_doc
        docs[class_name+'[all-methods]'] = class_doc + '\n' + self.create_method_markdown(False)
        docs[class_name+'[methods-with-docstrings]'] = class_doc + '\n' + self.create_method_markdown(True)

        return docs


    def create_method_markdown(self, only_with_docstrings=True):
        method_docs = []
        if self.methods:
            # reorder self.methods so __init__ comes first, if it isn't already
            methods = self.methods[:]
            for idx, meth in enumerate(self.methods):
                if meth.name == '__init__':
                    if idx != 0:
                        methods.remove(meth)
                        methods.insert(0, meth)
                    break

            for item in methods:
                if only_with_docstrings and not item.docstring:
                    continue
                if not check_excluded(item.name) and \
                   not check_excluded('{}.{}'.format(self.name, item.name)):
                        doc = item.create_markdown()
                        method_docs.append(doc)

        # join the methods into a single string
        method_docs = '\n'.join(method_docs)
        return method_docs


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
            if isinstance(doc, str):
                docs[item.name] = doc
            elif isinstance(doc, dict):
                docs.update(doc)
            else:
                raise RuntimeError("unknown type for doc")
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