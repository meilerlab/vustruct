# Identify files that lack exuecte by other permission
# A WORK IN PROGRESS
find . -type f '(' -user "mothcw" -perm -u=x ')' ! '(' -user "mothcw" -perm -o=x ')'

