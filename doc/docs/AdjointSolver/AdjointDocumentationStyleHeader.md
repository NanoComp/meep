<!----------------------------------------------------->
<!- LaTeX macros used only for adjoint documentation  ->
<!----------------------------------------------------->
$$
   \newcommand{\vb}{\mathbf}
   \newcommand{\wt}{\widetilde}
   \newcommand{\mc}{\mathcal}
   \newcommand{\bmc}[1]{\boldsymbol{\mathcal{#1}}}
   \newcommand{\sup}[1]{^{\text{#1}}}
   \newcommand{\pard}[2]{\frac{\partial #1}{\partial #2}}
   \newcommand{\VMV}[3]{ \Big\langle #1 \Big| #2 \Big| #3 \Big\rangle}
$$

<!----------------------------------------------------->
<!- custom CSS overrides for a singlee markdown page  ->
<!----------------------------------------------------->
<style type='text/css'>
  code,.rst-content { border:0;
                      font-size: 90%;
                    }

.highlight { background-color: #fdf6e3; color: #586e75 }
.highlight .c { color: #93a1a1 } /* Comment */
.highlight .err { color: #586e75 } /* Error */
.highlight .g { color: #586e75 } /* Generic */
.highlight .k { color: #859900 } /* Keyword */
.highlight .l { color: #586e75 } /* Literal */
.highlight .n { color: #586e75 } /* Name */
.highlight .o { color: #859900 } /* Operator */
.highlight .x { color: #cb4b16 } /* Other */
.highlight .p { color: #586e75 } /* Punctuation */
.highlight .cm { color: #93a1a1 } /* Comment.Multiline */
.highlight .cp { color: #859900 } /* Comment.Preproc */
.highlight .c1 { color: #93a1a1 } /* Comment.Single */
.highlight .cs { color: #859900 } /* Comment.Special */
.highlight .gd { color: #2aa198 } /* Generic.Deleted */
.highlight .ge { color: #586e75; font-style: italic } /* Generic.Emph */
.highlight .gr { color: #dc322f } /* Generic.Error */
.highlight .gh { color: #cb4b16 } /* Generic.Heading */
.highlight .gi { color: #859900 } /* Generic.Inserted */
.highlight .go { color: #586e75 } /* Generic.Output */
.highlight .gp { color: #586e75 } /* Generic.Prompt */
.highlight .gs { color: #586e75; font-weight: bold } /* Generic.Strong */
.highlight .gu { color: #cb4b16 } /* Generic.Subheading */
.highlight .gt { color: #586e75 } /* Generic.Traceback */
.highlight .kc { color: #cb4b16 } /* Keyword.Constant */
.highlight .kd { color: #268bd2 } /* Keyword.Declaration */
.highlight .kn { color: #859900 } /* Keyword.Namespace */
.highlight .kp { color: #859900 } /* Keyword.Pseudo */
.highlight .kr { color: #268bd2 } /* Keyword.Reserved */
.highlight .kt { color: #dc322f } /* Keyword.Type */
.highlight .ld { color: #586e75 } /* Literal.Date */
.highlight .m { color: #2aa198 } /* Literal.Number */
.highlight .s { color: #2aa198 } /* Literal.String */
.highlight .na { color: #586e75 } /* Name.Attribute */
.highlight .nb { color: #B58900 } /* Name.Builtin */
.highlight .nc { color: #268bd2 } /* Name.Class */
.highlight .no { color: #cb4b16 } /* Name.Constant */
.highlight .nd { color: #268bd2 } /* Name.Decorator */
.highlight .ni { color: #cb4b16 } /* Name.Entity */
.highlight .ne { color: #cb4b16 } /* Name.Exception */
.highlight .nf { color: #268bd2 } /* Name.Function */
.highlight .nl { color: #586e75 } /* Name.Label */
.highlight .nn { color: #586e75 } /* Name.Namespace */
.highlight .nx { color: #586e75 } /* Name.Other */
.highlight .py { color: #586e75 } /* Name.Property */
.highlight .nt { color: #268bd2 } /* Name.Tag */
.highlight .nv { color: #268bd2 } /* Name.Variable */
.highlight .ow { color: #859900 } /* Operator.Word */
.highlight .w { color: #586e75 } /* Text.Whitespace */
.highlight .mf { color: #2aa198 } /* Literal.Number.Float */
.highlight .mh { color: #2aa198 } /* Literal.Number.Hex */
.highlight .mi { color: #2aa198 } /* Literal.Number.Integer */
.highlight .mo { color: #2aa198 } /* Literal.Number.Oct */
.highlight .sb { color: #93a1a1 } /* Literal.String.Backtick */
.highlight .sc { color: #2aa198 } /* Literal.String.Char */
.highlight .sd { color: #586e75 } /* Literal.String.Doc */
.highlight .s2 { color: #2aa198 } /* Literal.String.Double */
.highlight .se { color: #cb4b16 } /* Literal.String.Escape */
.highlight .sh { color: #586e75 } /* Literal.String.Heredoc */
.highlight .si { color: #2aa198 } /* Literal.String.Interpol */
.highlight .sx { color: #2aa198 } /* Literal.String.Other */
.highlight .sr { color: #dc322f } /* Literal.String.Regex */
.highlight .s1 { color: #2aa198 } /* Literal.String.Single */
.highlight .ss { color: #2aa198 } /* Literal.String.Symbol */
.highlight .bp { color: #268bd2 } /* Name.Builtin.Pseudo */
.highlight .vc { color: #268bd2 } /* Name.Variable.Class */
.highlight .vg { color: #268bd2 } /* Name.Variable.Global */
.highlight .vi { color: #268bd2 } /* Name.Variable.Instance */
.highlight .il { color: #2aa198 } /* Literal.Number.Integer.Long */

.tomorrow-comment, pre .comment, pre .title {
  color: #8e908c;
}

.tomorrow-red, pre .variable, pre .attribute, pre .tag, pre .regexp, pre .ruby .constant, pre .xml .tag .title, pre .xml .pi, pre .xml .doctype, pre .html .doctype, pre .css .id, pre .css .class, pre .css .pseudo {
  color: #c82829;
}

.tomorrow-orange, pre .number, pre .preprocessor, pre .built_in, pre .literal, pre .params, pre .constant {
  color: #f5871f;
}

.tomorrow-yellow, pre .class, pre .ruby .class .title, pre .css .rules .attribute {
  color: #eab700;
}

.tomorrow-green, pre .string, pre .value, pre .inheritance, pre .header, pre .ruby .symbol, pre .xml .cdata {
  color: #718c00;
}

.tomorrow-aqua, pre .css .hexcolor {
  color: #3e999f;
}

.tomorrow-blue, pre .function, pre .python .decorator, pre .python .title, pre .ruby .function .title, pre .ruby .title .keyword, pre .perl .sub, pre .javascript .title, pre .coffeescript .title {
  color: #4271ae;
}

.tomorrow-purple, pre .keyword, pre .javascript .function {
  color: #8959a8;
}

pre code {
  display: block;
  background: white;
  color: #4d4d4c;
  font-family: Menlo, Monaco, Consolas, monospace;
  line-height: 1.5;
  border: 1px solid #ccc;
  padding: 10px;
}


</style>
