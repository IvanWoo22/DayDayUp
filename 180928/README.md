# CapQuiz
You have m dollars, and a bottle of Coca is l dollars. Now you can take n caps to exchange a new and full bottle.

How many bottles of Coca can you drink?


It just a little Perl program to kill time.






# Bulk Download sra Files

Just take [PRJEB8035](https://www.ebi.ac.uk/ena/data/view/PRJEB8035) as an example.

We can get a "TEXT" button at the bottom of this page. 
Affter click it, we can get a txt file named after `PRJEB8035.txt` which contains nearly all of this project's data infomation.

Then we can run like this:
```
Rscript BDsra.R ~/PRJEB8035.txt ~/PRJEB8035/
```
Here `~/PRJEB8035.txt` means the input file, while `~/PRJEB8035/` means the folder you want to download to.

The download utility we use in this script is [aria2](https://aria2.github.io/). It's lightweight multi-protocol & multi-source.
