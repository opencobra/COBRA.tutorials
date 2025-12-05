COBRA Tutorials
===============

Overview
---------
The COBRA.tutorials repository exists to collect and maintain tutorials for developed and published constraint-based genome-scale modelling methods and studies. Its purpose is to provide clear, reproducible instructional material that demonstrates how new algorithms, workflows, and analysis techniques can be applied within the COBRA modelling framework. By centralising these tutorials, the repository supports both new and experienced users in learning, comparing, and adopting state-of-the-art constraint-based modelling approaches.

Repository structure
---------------------
Tutorials are here to get you started with using `The COBRA
Toolbox <https://opencobra.github.io/cobratoolbox>`__. The
tutorials are grouped according to the ``src/`` folder structure:

- `analysis <https://github.com/opencobra/COBRA.tutorials/tree/master/analysis>`__
- `base <https://github.com/opencobra/COBRA.tutorials/tree/master/base>`__
- `dataIntegration <https://github.com/opencobra/COBRA.tutorials/tree/master/dataIntegration>`__
- `design <https://github.com/opencobra/COBRA.tutorials/tree/master/design>`__
- `reconstruction <https://github.com/opencobra/COBRA.tutorials/tree/master/reconstruction>`__
- `visualization <https://github.com/opencobra/COBRA.tutorials/tree/master/visualization>`__

All tutorials are provided in these formats: ``.mlx``, ``.m``, and ``.html``.

- The interactive version ``.mlx`` is a MATLAB Live Script format and can be run using `the MATLAB Live-script editor <https://nl.mathworks.com/help/matlab/matlab_prog/what-is-a-live-script.html>`__.
- The static version ``.html`` is automatically generated and is accessible on the `tutorial section of the COBRA Toolbox documentation <https://opencobra.github.io/COBRA.tutorials>`__.
- The ``.pdf`` version can be downloaded from the same tutorial section. The ``.m`` version can be opened and run directly in MATLAB, which is particularly useful for building new analysis scripts based on existing tutorials.

How the Continuous Integration (CI) System Works
------------------------------------------------

The COBRA.tutorials repository uses an automated **continuous integration (CI)** workflow.  
When a contributor pushes a new or updated ``.mlx`` tutorial to the repository:

1. The CI pipeline is automatically triggered.  
2. The ``.mlx`` file is converted into three formats:  
   - ``.m`` (MATLAB script)  
   - ``.pdf`` (print-ready copy)  
   - ``.html`` (web-friendly version)  
3. These generated files are published to the COBRA Toolbox website, making the tutorial immediately accessible to users.  

The diagram below illustrates this workflow:

.. raw:: html

   <div style="text-align: center; width: 100%; margin-top: 10px; margin-bottom: 20px;">
     <img src="doc/img/COBRA_Tutorials_CI_Workflow.png"
          alt="COBRA Tutorials CI Workflow"
          style="display: block; margin-left: auto; margin-right: auto; width: 80%; max-width: 900px;">
   </div>

Contribute a new tutorial or modify an existing tutorial
========================================================

A template for generating a new tutorial is provided `here
<https://github.com/opencobra/COBRA.tutorials/tree/master/.template/tutorial_template.mlx>`__.

There are two ways to contribute to the tutorials:

A) Contribute using ``git`` (via command line)
----------------------------------------------

Fork and checkout your branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Fork the `COBRA.tutorials repository <https://www.github.com/opencobra/COBRA.tutorials>`__ on Github.

2. Clone the forked repository to a directory of your choice:

   .. code-block:: console

      git clone git@github.com:<userName>/COBRA.tutorials.git fork-COBRA.tutorials.git

3. Change to the directory:

   .. code-block:: console

      cd fork-COBRA.tutorials.git/

4. Set the upstream to the ``opencobra/COBRA.tutorials`` repository:

   .. code-block:: console

      git remote add upstream git@github.com:opencobra/COBRA.tutorials.git

5. Fetch from the upstream repository:

   .. code-block:: console

      git fetch upstream

6. Checkout a new branch from ``upstream/master``:

   .. code-block:: console

      git checkout -b master upstream/master

7. Now, make your changes to the tutorial in MATLAB.

Submit your changes and open a pull request to the ``master`` branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

8. Once you are done making changes, add the files to your branch (``tutorial_<yourFile>`` refers to the tutorial name). Make sure to add the ``.mlx`` format of the tutorial.

   .. code-block:: console

      git add tutorial_<yourFile>.mlx
      git commit -m "Changes to tutorial_<yourFile>"

9. Push your commits on ``<yourBranch>`` to your fork:

   .. code-block:: console

      git push origin <yourBranch>

10. Browse to your fork on:

    * ``https://www.github.com/<yourUserName>/COBRA.tutorials``

11. Click on **Compare & Pull Request**.

12. Confirm the target branch is ``master``.

13. Submit your pull request.

14. Wait until your pull request is accepted.

.. |icon_analysis| raw:: html

   <img src="https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable/_static/img/analysis.png" height="14px">

.. |icon_base| raw:: html

   <img src="https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable/_static/img/base.png" height="14px">

.. |icon_dataIntegration| raw:: html

   <img src="https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable/_static/img/dataIntegration.png" height="14px">

.. |icon_design| raw:: html

   <img src="https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable/_static/img/design.png" height="14px">

.. |icon_reconstruction| raw:: html

   <img src="https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable/_static/img/reconstruction.png" height="14px">

.. |icon_visualization| raw:: html

   <img src="https://github.com/opencobra/cobratoolbox/tree/gh-pages/stable/_static/img/visualization.png" height="14px">

B) Contribute using the ``MATLAB.devTools``
-------------------------------------------

You can use the `MATLAB.devTools <https://github.com/opencobra/MATLAB.devTools>`__ to submit your tutorial.
