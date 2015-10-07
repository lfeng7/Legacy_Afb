read  -e -p "Enter the commit message: " commit_message
read -e -p "Enter the branch name for push: " branch_name
git add *.*
git commit -m "'$commit_message'"
git push origin "$branch_name"
