.. scGEAToolbox documentation master file, created by
   sphinx-quickstart on Sun Dec 20 17:44:55 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scGEAToolbox's documentation!
========================================

scGEAToolbox is a Matlab toolbox for scRNA-seq data analysis.

#### Official Websites and Social Networks

Please, visit the official website [![View scGEAToolbox on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox) of scGEAToolbox for further information. 

You might also want to take a look at scGEAToolbox in action: see the [Youtube playlist](https://youtube.com/playlist?list=PLekWsTmNpdrHNznPesjE7dWxx7iMrucTo)

Visit the dedicated [Slack Channel](https://join.slack.com/t/scgeatoolbox/shared_invite/zt-6zl1893a-S4WQKH3XPb_q68r0ejoMEA) if you have questions or to report bugs.

Instead, if you create amazing visualizations using scGEAToolbox and you want to tweet them, remember to include the "#scGEAToolbox" hashtag: we will be happy to retweet.

You might be interested in taking a look at the [official documentation of scGEAToolbox](https://scgeatoolbox.readthedocs.io/), with posts about applications (analysis, visualization or both) to cool datasets.

Getting Started
===============

* Easy ability to install/use as a Matlab toolbox (tip o' the hat to 
  `James Cai's scGEAToolbox <https://github.com/jamesjcai/scGEAToolbox>`_ for
  showing the way);
* Style tweaks compared to the source themes, such as better code-block
  alignment, Github button placement, page source link moved to footer,
  improved (optional) related-items sidebar item, and many more;
* Many customization hooks, including toggle of various sidebar & footer
  components; header/link/etc color control; etc;
* Improved documentation for all customizations (pre-existing & new).

Quick demos
-----------

Read the Docs offers many advanced features and options.
Learn more about these integrations and how you can get the most
out of your documentation and Read the Docs.

First steps
-----------

Are you new to software documentation
or are you looking to use your existing docs with Read the Docs?
Learn about documentation authoring tools such as Sphinx and MkDocs
to help you create fantastic documentation for your project.

Step-by-step Guides
-------------------

These guides will help walk you through specific use cases
related to Read the Docs itself, documentation tools like Sphinx and MkDocs
and how to write successful documentation.

* **Weighted List**:
  :doc:`With Sphinx </intro/getting-started-with-sphinx>` |
  :doc:`With MkDocs </intro/getting-started-with-mkdocs>` |
  :doc:`Feature Overview </features>` |
  :doc:`/choosing-a-site`
  
Project background
==================

scGEAToolbox is a modified (with permission) version of `Kenneth Reitz's
<http://kennethreitz.org>`_ `"krTheme" Sphinx theme 
<https://github.com/jamesjcai/scGEAToolbox>`_ (it's the one used 
in his `Requests <http://python-requests.org>`_ project). Kenneth's 
theme was itself originally based on Armin Ronacher's `Flask
<http://flask.pocoo.org/>`_ theme. Many thanks to both for their hard work.


Implementation notes
====================

* `Fabric #419 <https://github.com/fabric/fabric/issues/419>`_ contains a lot of
  general exposition & thoughts as I developed scGEAToolbox, specifically with a
  mind towards using it on two nearly identical 'sister' sites (single-version
  www 'info' site & versioned API docs site).
* scGEAToolbox includes/requires a tiny Sphinx extension on top of the theme
  itself; this is just so we can inject dynamic metadata (like scGEAToolbox's own
  version number) into template contexts. It doesn't add any additional
  directives or the like, at least not yet.

Contents:
=========

.. toctree::
   :maxdepth: 2

   getting_started
   case_studies
   case_studies2
   readme

