This folder contains various latex-templates. They all need data from the post-processing produced by calling rungeneric.py on experimental data (see http://coco.lri.fr/COCOdoc/firsttime.html). All templates use as default input folder ppdata (post-processed "data"). Which template to use depends on the testbed that has been run to do the experiments and how many algorithms have been used with rungeneric.py. 

Usage example:

   pdflatex templateBBOBarticle.tex


Templates for GECCO BBOB:

noise-free testbed: 
  templateBBOBarticle.tex - for one algorithm (rungeneric.py called with one argument)
  templateBBOBcmp.tex - for two algorithms (rungeneric.py called with two arguments)
  templateBBOBmany.tex - for more than two algorithms (more than two arguments calling rungeneric.py)
  
noisy testbed:
  templateBBOBnoisyarticle.tex - for one algorithm (rungeneric.py called with one argument)
  templateBBOBnoisycmp.tex - for two algorithms (rungeneric.py called with two arguments)
  templateBBOBnoisymany.tex - for more than two algorithms (more than two arguments calling rungeneric.py)

Note: The same template can be used if run with low ("expensive") or "normal" budget (figures, tables and captions are adapted automatically during the postprocessing)

Non-specific templates:

template1generic.tex - generic template for one algorithm on the noise-free testbed
template2generic.tex - generic template for two algorithms on the noise-free testbed
template3generic.tex - generic template for more than two algorithms on the noise-free testbed

noisytemplate1generic.tex - generic template for one algorithm on the noisy testbed
noisytemplate2generic.tex - generic template for two algorithms on the noisy testbed
noisytemplate2generic.tex - generic template for more than two algorithms on the noisy testbed




