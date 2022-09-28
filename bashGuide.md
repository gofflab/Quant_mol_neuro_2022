# Bash Cheat sheet

Find where you are: `pwd`
What is in here: `ls`

- A program is an _executable_ file.
- Any command you call in bash is a program and it has to be somewhere in the system. Where bash looks is called `PATH`.
To find where that program is, use `which`.
- If you have a program in your folder, you must prepend `./` to the program so that shell knows that the program is here.

https://www.rozmichelle.com/pipes-forks-dups/
- Programs take in something and spit out something (although they might do something else in the background).
- Two channels that programs take in.
	- Arguments
		- These are all the words we put in behind the program. Usually specify some configurations. Ordering matters.
		- Arguments without dashes are usually required.
		- An argument with a double dash is usually an option and may/may not have extra info. Equal signs can be use or not. Ex. `--color=red` or `--color red`.
    Another one is `--fast`. This one doesn't have extra info.
		- An argument with a dash is always a single character and usually an abbreviation of the double dash. Ex. `-c red`.
    If you see a single dash with multiple characters, ex. `-cd`, that usually means `-c -d`.
	- STDIN
		- A channel that data flows in.

- Two channels where data flow out.
	- STDOUT
	- STDERR
		- These are effectively the same but separated for organizational reasons. By default, these just flow into the terminal.
		- STDOUT is usually the data.
		- STDERR is usually the log or notes the program made.

		- To save STDOUT to a file, use >. Ex. `ls > file.out` saves the output of ls to a file named file.out.
		- To save STDERR to a file, use 2>. Ex. `ls 2> file.out` creates file.out, which has nothing inside because ls doesn't output STDERR.

- With this, you can chain programs. For example, 
	`ls | head -n 1`
  `ls` outputs what we have in our current directory through STDOUT, which is passed into STDIN via the pipe '|' operator of `head`.
  Using the argument -n 1, `head` spits out through STDOUT only the first line it received.

Go here: cd
	- Absolute paths always start with `/`.
	- Your home directory (think of it as desktop) is at `/home/$USERNAME`
	- Relative paths (not started with `/`) are related to where you are, use `pwd`.
	- `.` (the dot) means current directory.
	- `..` (two dots) mean the parent directory (directory above).
	- To go to a sibling folder, you can do `cd ../folder/`
	- Directory paths can end with `/` but not files.

I don't want to type this again and again. Use variables. = defines a variable.
You don't need to quote the thing you're setting if you're not using spaces or some weird characters.
For example,
	`VAR=Hello` is fine.
	`VAR=Hello world` is not. Use `VAR="Hello world"`.

How do we discern literal texts and variables? The $ dollar sign.
Example: suppose we set F=folder
	`cd folder` is the same as `cd $F`.

Strings can be concatenated or stick together using braces {}.
- `${F}${F}${F}` => folderfolderfolder
- `${F}here` => folderhere
- `$Fhere` => Nothing. Bash tried to look for a variable called `Fhere`.
