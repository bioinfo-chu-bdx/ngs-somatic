AGeNT
=======
The <strong>A</strong>gilent <strong>Ge</strong>nomics <strong>T</strong>oolKit (AGeNT) is a collection of command-line tools, developed by Agilent, for NGS data processing. With minimal configuration, you can start using these tools on the Linux or Windows command-line.
<br>
## Command-line syntax

### Getting help with AGeNT
To get help with any command while using AGeNT, simply type help at the end of any command name.

For example, this command displays help for the general AGeNT options and the available top-level commands:
<code>$ agent help</code>

The following command displays the specific help message for the LocatIt module:
<code>$ agent locatit help</code>

### Command Structure in AGeNT

AGeNT uses a multipart structure on the command line that must be specified in this order:
1. The base call to the AGeNT program.
2. The top-level command, which typically corresponds to a tool supported by AGeNT.
3. The subcommand that specifies which operation to perform.
4. General CLI options or parameters required by the operation.

Structure:
<code>$ agent \<command\> \<subcommand\> [subcommand-specific options and parameters] 
</code>

Input parameters can consist of any arguments supported by the AGeNT modules, such as numbers, strings, lists, maps, and JSON structures. The supported arguments depend upon the specific command and subcommand specified.

### Example
<code>
$ agent locatit -v2Only -d 1 -m 3 -q 25 \<br>
-l /Users/uname/data/Covered.bed \<br>
-o /Users/uname/data/test_output.bam \<br>
/Users/uname/data/test_input.bam
</code>

<br>

**Note:** *Paths and filenames should also contain no spaces or other non-permissible characters on a Unix or Windows command-line.*

