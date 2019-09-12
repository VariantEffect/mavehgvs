# hgvs-patterns
A python utility containing HGVS Regex patterns to match a subset of the 
[HGVS](http://varnomen.hgvs.org/) standard. 

# Supported Syntax
Regex patterns have been implmeneted for DNA, RNA and Protein **subsitutions**, **deletions**, **insertions**, 
**deletion-insertions** and **frame-shifts**. Support for all syntax in the HGVS guidelines
within these events is not currently implemented. In general, patterns matching 
the [examples](http://varnomen.hgvs.org/recommendations/DNA/variant/substitution/) given on each
event page are matchable.

# Installation
To install `hgvs-patterns` use either of the commands below

```bash
pip install git+https://github.com/VariantEffect/hgvs-patterns.git
```

Or

```bash
git clone https://github.com/VariantEffect/hgvs-patterns
cd hgvs-patterns
python setup.py install
```

# Getting Started
To use hgvs-patterns in your project you can import regex patterns from
the top level to match any type of variant.

```
>>> from hgvsp import single_variant_re, multi_variant_re

>>> single_variant_re.fullmatch('p.Arg78_Gly79ins23')
<_sre.SRE_Match object; span=(0, 18), match='p.Arg78_Gly79ins23'>

>>> multi_variant_re.fullmatch('c.[123A>G;125G>T]')
<_sre.SRE_Match object; span=(0, 17), match='c.[123A>G;125G>T]'>
``` 

Alternatively, you can match specific events using the `dna`, `rna` and `protein`
modules. Each of the modules `dna`, `rna` and `protein` follow the same import patterns
and define regex expressions using the same variable names. Patterns containing
`_variant_` in their name will only match mutation events prefixed with a HGVS prefix such
as those from `cnpmgr`. Those with `_event_` in their name will only match specific
events without prefixes. All other regex patterns will match events with or without
prefixes.

```
>>> from hgvsp.dna import deletion_re, insertion_re, single_variant_re, any_event_re

>>> single_variant_re.fullmatch('123_124delinsAAA')
None
>>> single_variant_re.fullmatch('n.123_124delinsAAA')
<_sre.SRE_Match object; span=(0, 18), match='n.123_124delinsAAA'>

>>> deletion_re.fullmatch('c.123del')
<_sre.SRE_Match object; span=(0, 17), match='c.123del'>

>>> insertion_re.fullmatch('123_124insGATTACA')
<_sre.SRE_Match object; span=(0, 17), match='123_124insGATTACA'>

>>> any_event_re.fullmatch('c.123_124insGATTACA')
None
>>> any_event_re.fullmatch('123_124insGATTACA')
<_sre.SRE_Match object; span=(0, 17), match='123_124insGATTACA'>
``` 

