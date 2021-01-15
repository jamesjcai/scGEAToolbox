.. _getting_started:

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

scGEAToolbox is a MATLAB-based tool for analyzing and visualizing scRNA-seq data as interactive high-dimensional exploratory data analysis tool (see :ref:`intro_heatmap_clustergram`). This section will provide information on how to interact with the visualization and how to quickly visualize your own data using the :ref:`scGEAToolbox_web`, :ref:`scGEAToolbox2`, and :ref:`scGEAToolbox_widget`.

See :ref:`case_studies` for examples of how scGEAToolbox can be used to explore and analyze real world data. For developers interested in building their own functions using scGEAToolbox, please refer to the :ref:`building_web_page` section.

What's New
=============

scGEAToolbox YouTube Channel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">


         <iframe width="560" height="315" src="https://www.youtube.com/embed/qHb-RtvN4D4" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

We recently presented.

