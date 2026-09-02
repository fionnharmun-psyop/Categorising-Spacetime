When you have downloaded Julia you want to download these files to a folder then open that folder as a new environment in VS code. Having done so, you then open the Julia REPL (it's like a terminal).
It should be green and have julia>. Then press ] in the REPL to enter package mode. It should go blue with (@v1.12) pkg> meaning you are in the global Julia 1.12 environment on your PC. Enter: activate "." exactly.
This tells Julia to focus on your current environment, and the REPL icon should change to (YourFolder) pkg>. Then type instantiate. This tells Julia to look and the manifest and project files and to download them from
the internet. The manifest is preloaded with the Plots and Gradus and other Packages. Doing it this way stops any version dependency conflicts. Say, if Gradus needs Plots V2 to work but the latest version is Plots V3, 
Julia will detect that and load Plots V2 so it all works. It will do a lot of downloading and pre-compiling then after you can run the example Gradus script. To import more packages make sure you are in the working
environment, then, in package mode, enter add YourPackage. After, you must push these changes onto GitHub with the script and new manifest and project files so that when another person pulls the changes everything works well.

