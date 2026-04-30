#### Martín Santamarina García
#### 05/08/2025

######################
#### GIT commands ####
######################

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
