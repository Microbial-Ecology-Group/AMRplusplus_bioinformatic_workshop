## Training portal

To begin the hands-on portion of the workshop we will load up the training portal. We recommend using Google Chrome due to copy and paste feature that works well with the portal.

Paste this web address into your web browser:

```

portal-lms.hprc.tamu.edu

```

As seen in the image below the dropdown box is where you will place your username and password.


![Portal login](../resources/portal-images/portal_login.png)


Now that you are logged into your portal dashboard...Click on the `Files` dropdown tab and you will be redirected into your training `home` directory.


![Portal dashboard](./resources/portal-images//portal_fileTab.png)


Here, you will experience a GUI file explorer, this is where we can use the feature to `Upload` or `Download` files to the portal. This will come in handy to download and explore locally, the count matrix output files produced by AMR++.


![Portal files](./resources/portal-images/portal_dashHomeDir.png)


If you are in need of support from our "behind-the-scenes" instructors, we will navigate to the `My Interactive Sessions` tab (This will become avialable once we load up a Interactive Desktop in the Bioinformatics section). Once on this page we will right click on the `View only (Share-able link)`.


![Portal view only](./resources/portal-images/portal_viewOnly.png)


This will pull up a new Interactive Desktop tab (as seen in the image below). You will then copy the URL from your web browser and paste this into the ASM Slack group to either Dr. Perez (Lisa Perez) or Dr. Pinnell (ljpinnell). They will now be able to have a live view your Interactive Desktop and be able to help you troubleshoot any problems going on!


![Portal view only tab](./resources/portal-images/portal_viewOnlyLink.png)


When you are ready to logout of the portal, click the `Log Out` tab (top right of your browser). Once loaded, you will then see a message on how to completely log out of the training portal.


![Portal logout](./resources/portal-images/portal_logout.png)


## Bioinformatics

Now from your dashboard we will launch the interactive compute nodes to connect us to the HPC cluster resources. This is where we will run the AMR++ pipeline! Hover over the `Interactive Apps` dropdown tab and click on the `mm Desktop` button to connect.


![Portal interactive tab](./resources/portal-images/portal_appTab.png)


You will now be able to choose the number of hours and cores you want to use for the workshop. Click the `Launch` to start the queue to load the Interactive Desktop.


![Portal resources](./resources/portal-images/portal_clusterRescs.png)


You will be queued for a brief moment and then once your session is created we will click the `Launch mm Desktop` button to launch a new web browser tab with your Interactive Desktop.


![Portal node](./resources/portal-images/portal_node.png)


Now on your new Interacivte Desktop tab, click on the terminal emulator app (the second black app at the bottom of the desktop) and then we will begin the command-line portion of the workshop!


![Portal terminal](./resources/portal-images/portal_termEm.png)


The AMR++ pipeline, the important software dependencies, and configuration files are are pre-installed and ready for you to run on the portal! We will need to do some basic navigation and listing of the file structure before we run the AMR++ pipeline. It is important to note the `$` character in the code block examples below represents the end of the command prompt, therefore you will not need to type in or copy it when you run your commands.

The first command we will run is `pwd`, this stands for "print working directory" and it will print out the `absolute path` of where you are currently standing on the file system.
```bash
$ pwd
```

To change directories we will use the `cd` command and we will give it again the `aboslute path` with the env variable `$HOME` to your `Desktop` directory.
```bash
$ cd $HOME/Desktop/
```

The `ls` command without any options or arguments will output the files and directories of where you are currently standing (i.e. your Desktop).
```bash
$ ls
```

Here you will now see all the files that you will need to run AMR++ and files needed to run the hands-on statistics part of the workshop. The files needed for the AMR++ run will be the bash script `run_AmrPlusPlus.sh` and the nextflow congfiguation file `nextflow.config`.


![Portal commands](./resources/portal-images/portal_basicCmds.png)


To start your AMR++ run we will run the following command:


```bash
$ bash run_AmrPlusPlus.sh
```

The standard output to the terminal screen is interactive. Once the run has finished it should look similar as the image below:


![Portal AMR out](./resources/portal-images/portal_AMRout.png)


The directory with the AMR++ results will be outputted here in your `Desktop` directory where you are currently standing. It will be named `Small_Results_WithKraken0.1`

There will be several sub-directories within our parent output directory. Within the sub-directories `ResistomeResults` and `KrakenResults` we will locate our MEGARes resistome and kraken2 microbiome count matrix CSV files respectively. There are also other important output stat files from the run e.g. the `trimmomatic.stats` file and the `host.removal.stats` file. The image below is a condensed version of the AMR++ output directory structure.


![Portal output tree](./resources/portal-images/portal_amrOutTree.png)


Next we will move back to our portal dashboard and download our AMR++ output files locally to our computers to further explore them!


![Portal download](./resources/portal-images/portal_dwnl.png)


## Statistics

To start the statistics portion of the workshop we will connect back to our Interactive Desktop node. Like before we will navigate to our `Desktop` directory. Here is where we will launch the Rstudio application.

```bash
$ rstudio
```
On the bottom right panel of Rstudio in the `Files` tab click on the file `AMR_stats_portal.R`. This will pull up the R script we will use to analyze our metagenomic output files. There are a few different ways to step through the code. You can use the keystroke `CTRL+Enter` or `CTRL+Return`. Or you can use the `Run` button at the top of the script. You can also highlight or select specific line(s) you want to run with the keystroke or run button.


![Portal Rstudio](./resources/portal-images/portal_r.png)


While running the running the R code, if you see the warning below, no worries this may just be due to your computer screen dimensions and the portal.


```
Warning message:
In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  X11 used font size 8 when 9 was requested
```
