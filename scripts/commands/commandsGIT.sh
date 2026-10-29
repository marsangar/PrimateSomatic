#### Martín Santamarina García
#### 05/08/2025

######################
#### GIT commands ####
######################


1. One-time setup: set your identity

git config --global user.name "marsangar"
git config --global user.email "marsangar15@gmail.com"

2. Create a new repo locally
mkdir my-project && cd my-project
git init -b main

3. Add content and make your first commit
echo "# My Project" > README.md
git add .
git commit -m "Initial commit"

4. Connect to a remote and push
git remote add origin git@github.com/marsangar/xsmut.git
git remote set-url origin git@github.com:marsangar/xsmut.git
git push -u origin main

5. Everyday cycle: update, commit, push
git status
git diff
git add path/to/file.txt
git commit -m "Describe what changed"
git pull --rebase
git push

Useful extras
Task	Command
View history	git log --oneline --graph
Create and switch to a branch	git switch -c feature-x
Push a new branch	git push -u origin feature-x
Unstage a file	git restore --staged file.txt
Discard local edits to a file	git restore file.txt
See remotes	git remote -v
Clone an existing repo	git clone <url>

Tips:

Use the HTTPS URL (https://github.com/USERNAME/my-project.git) if you haven't set up SSH keys. GitHub then asks for a personal access token, not your password.
Add a .gitignore file to keep files out of the repo, such as .DS_Store, __pycache__/ and large data files.








### This script is used to keep track of relevant git commands

git status
git add .
git commit -m "Your commit message here"
git push origin main




#### To generate a new SSH key pair, you can use the following command:
#Step 1 — generate SSH key
ssh-keygen -t ed25519 -C "marsangar15@gmail.com"

#Step 2 — add key to ssh-agent
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519

 #Step 3 — copy public key
cat ~/.ssh/id_ed25519.pub

#Step 4 — change repo remote URL
#Check current remote:
git remote -v
#Switch to SSH:
git remote set-url origin git@github.com:marsangar/PrimateSomatic.git

#Step 5 — test connection
ssh -T git@github.com



####
Why SSH is better for your project

For PrimateSomatic (HPC + multi-machine workflow):

✔ no tokens to expire
✔ works on HPC clusters easily
✔ no repeated login prompts
✔ more stable for pipelines and automation
✔ standard in bioinformatics labs
