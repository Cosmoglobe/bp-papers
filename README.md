## BeyondPlanck Papers

The full list of BP Papers (in PDF format) appear at:

https://papers.beyondplanck.science

Since these are unreleased papers, they are not publicly available and are behind a:

`bps3cr3t.` (note the period at the end)

**NOT** all PDFs are properly compiled yet (read on)

There is a Gitlab process that compiles all PDFs, once a commit to the repo takes place.
This is done "automagically", nothing is required from your part.

## General Instructions

How to make sure your paper is properly compiled:

- There is an execution pipeline that displays the current compilation tasks at:
  https://gitlab.com/BeyondPlanck/papers/pipelines (also accessible from CI/CD ->
  Pipelines menu)

- It displays a Status and the User that triggered it (more likely you)

- The status is clickable and it shows you all the executing Jobs.

- The Jobs are split in multiple 'paper_XX' jobs which are the latex compilation task
  for each paper (numbered as in HKE's overview paper)

- Each Job is appropriately colored (Green: ok, orange/red: faulty) and they are
  clickable.

- Clicking an a Job (more likely of your paper) shows you the latex compilation log, and
  potentially a useful error on why it failed.

- If your job fails to compile correctly, it will not be available in the
  https://papers.beyondplanck.science page (obviously)

- Fix the error and recommit (read Notes for help)

- The whole pipeline takes around 4 minutes to complete (as of now, with NOT all papers
  compiling)

## Notes

- The papers are compiled with 'latexmk'. Actual command line is:

`latexmk -cd -pdf -output-directory=../builds aDir/paper.tex -jobName=BP_XX`

executed from the root of the 'papers' repo. Parameters:

- `-cd`: Changes to the given dir while compiling
- `-pdf`: Generates PDF
- `-output-directory`: Outputs ALL papers in a single dir (builds/)
- `-jobName`: Names each compiled paper as BP_XX

* Paper authors are NOT using a standardized way of compiling their papers. I made the
  assumption that 'latexmk' ought to compile them all (and it does in the 5-6 I have
  tried so far).

* You should be able to run the same command locally on your machine to see if it
  compiles and if it causes the same error as on Gitlab.

* If you would rather use an alternate way of compiling (although we SHOULD be as common
  as possible) let me know the command line to compile your paper differently and I'll
  change the compilation process

## Proposal for Changes/Improvements

- Paper directories are hit and miss and very long (making script references difficult
  and very long). Can we all please, use:

  `paper_XX`

- Can we already use the final name of the PDF file? I do not know if there is a
  specific pattern that all Planck (or BeyondPlanck now) papers are using. Maybe
  "BP_Paper_title_in_camel_case"? As of now the names are all over the place, and I
  produce them all, temporarily for now, as `BP_XX.pdf`

- If you choose to help out and tidy things up (please DO), in order to do so without
  breaking the compilation pipeline, you will also have to change folder name and paper
  name in a couple of files too:

  1. `.gitlab-ci.yml`:

     - Change the relevant 'paper_xx' stanza with the new folder and paper name. (the
       'aDir/paper.tex' and the 'jobName' from the command above)

     - If your `paper_XX` stanza is not there (I run out of steam creating them all) you
       can copy and paste an existing one (editing accordingly) and create your own.

     - If you DO add your own, add a dependency to the 'website' section too (you'll see
       it)


    1. `website/index.html`:

        * You can edit your Paper Title as it displays in the homepage and the name of
          the compiled PDF there (it's plain HTML, it should be self explanatory) For
          href=pdf_name, no folder required, where pdf_name should be the same as the
          jobName from the latexmk command above

That should be all...

It might be too much info but all of this is actually boilerplate code and only need to
be set up once. Then it ought to work on it's own.

Do not hesitate to get in touch with me (@stratosgear gerakakis@planetek.gr) for any
further clarifications.

#### Footnote:

The total time required for the papers to be compiled and uploaded to the website
depends on the number of available Gitlab runners that do the processing. Gitlab offers
a limited amount of computing minutes per project per month, that eventual will run out
before each month's end. I have thus also added a runner that runs on a private personal
server of mine, that should always by available to take the remaining load.

If anyone has access to an always on server, capable of running Docker, and willing to
donate some spare cycles [do let me know](mailto:gerakakis@planetek.gr) for instructions
on how to join the Gitlab pool (hint: read
[gitlab-runner-readme.md](gitlab-runner-readme.md))

A Gitlab runner is also able to execute locally on a laptop, even for a few hours per
day, if your remember to start it. If even some of us were running one, that would be
more than enough to cover our needs.
