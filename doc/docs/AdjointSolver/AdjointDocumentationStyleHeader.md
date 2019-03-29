<!----------------------------------------------------->
<!- LaTeX macros used only for adjoint documentation  ->
<!----------------------------------------------------->
$$
   \newcommand{\vb}{\mathbf}
   \newcommand{\wt}{\widetilde}
   \newcommand{\mc}{\mathcal}
   \newcommand{\bmc}[1]{\boldsymbol{\mathcal{#1}}}
   \newcommand{\sup}[1]{^{\text{#1}}}
   \newcommand{\sups}[1]{^{\text{#1}}}
   \newcommand{\sub}[1]{_{\text{#1}}}
   \newcommand{\subs}[1]{_{\text{#1}}}
   \newcommand{\pard}[2]{\frac{\partial #1}{\partial #2}}
   \newcommand{\VMV}[3]{ \Big\langle #1 \Big| #2 \Big| #3 \Big\rangle}
$$

<style>
.superfences-tabs {
  display: flex;
  position: relative;
  flex-wrap: wrap;
}

.superfences-tabs .highlight {
  background: #ddd;
}

.superfences-tabs .superfences-content {
  display: none;
  order: 99;
  width: 100%;
}

.superfences-tabs label {
  width: auto;
  margin: 0 0.5em;
  padding: 0.25em;
  font-size: 120%;
  cursor: pointer;
}

.superfences-tabs input {
  position: absolute;
  opacity: 0;
}

.superfences-tabs input:nth-child(n+1) {
  color: #333333;
}

.superfences-tabs input:nth-child(n+1):checked + label {
    color: #FF5252;
}

.superfences-tabs input:nth-child(n+1):checked + label + .superfences-content {
    display: block;
}
</style>
