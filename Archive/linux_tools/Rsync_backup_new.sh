##Recover
#rsync --compress --archive --verbose --exclude={"dev","proc","sys","tmp","run","mnt","media","lost+found","cloud_research"} --hard-links --human-readable --inplace --numeric-ids --delete -e ssh /home/back_up/Tang/Daily_Backup/05252019/ root@tang.phys.gwu.edu:/ > /home/back_up/Tang/Tang.recover.$(date +%m%d%Y_%H.%M).log



#!/bin/sh
#
# This script is used on a QNAP TS-269 PRO. https://www.en0ch.se/qnap-and-rsync/
# 
# You have to change:
# 1. $SHAREUSR
# 2. $EXCLUDES (if you want o change the name of the file servername.excludes)
# 3. $SOURCE & $DESTINATION
# 4. user@yourserver.se for the mysqldump 
# 5. --password=SUPERSECRET
 
TODAY=$(date +"%Y%m%d")
LASTBACKUP=`date -d "7 day ago" +"%Y%m%d"`
LASTBACKUP="06062019"


# Set the path to rsync on the remote server so it runs with sudo.
#RSYNC="/usr/bin/sudo /usr/bin/rsync"

# Set the folderpath on the QNAP
# Dont't forget to mkdir $SHAREUSR
#SHAREUSR="/share/CACHEDEV1_DATA/yourserver.se"
 
# This is a list of files to ignore from backups.
# Dont't forget to touch $EXCLUDES
#EXCLUDES="$SHAREUSR/servername.excludes"

#LOG file
# Dont't forget to touch $LOG
LOG="/home/crontab_log/Tang.daily.backup"
 
# Remember that you will not be generating
# backups that are particularly large (other than the initial backup), but that
# you will be creating thousands of hardlinks on disk that will consume inodes.

#Folder To be backup and Destination

Folder_Backup="/home"

Destination_SSH="xli@128.164.54.240"

Destination_Path="/home/back_up/Weekly_Backup"

Destination="${Destination_SSH}:${Destination_Path}/${TODAY}"


 
# Keep database backups in a separate directory.
#mkdir -p $SHAREUSR/$TODAY/db

# This command rsync's files from the remote server to the local server.
# Flags:
#   -z enables gzip compression of the transport stream.
#   -e enables using ssh as the transport prototcol.
#   -a preserves all file attributes and permissions.
#   -x (or --one-file-system) Don’t cross filesystem boundaries
#   -v shows the progress.
#   --rsync-path lets us pass the remote rsync command through sudo.
#   --exclude-from points to our configuration of files and directories to skip.
#   --numeric-ids is needed if user ids don't match between the source and
#       destination servers.
#   --delete -r(ecursive) Deletes files from $DESTINATION that are not present on the $SOURCE
#   --link-dest is a key flag. It tells the local rsync process that if the
#       file on the server is identical to the file in ../$YESTERDAY, instead
#       of transferring it create a hard link. You can use the "stat" command
#       on a file to determine the number of hard links. Note that when
#       calculating disk space, du includes disk space used for the first
#       instance of a linked file it encounters. To properly determine the disk
#       space used of a given backup, include both the backup and it's previous
#       backup in your du command.
#
# The "rsync" user is a special user on the remote server that has permissions
# to run a specific rsync command. We limit it so that if the backup server is
# compromised it can't use rsync to overwrite remote files by setting a remote
# destination. I determined the sudo command to allow by running the backup
# with the rsync user granted permission to use any flags for rsync, and then
# copied the actual command run from ps auxww. With these options, under
# Ubuntu, the sudo line is:
#
#   rsync	ALL=(ALL) NOPASSWD: /usr/bin/rsync --server --sender -logDtprze.iLsf --numeric-ids . /
#
# Note the NOPASSWD option in theDestination_Path sudo configuration. For remote
# authentication use a password-less SSH key only allowed read permissions by
# the backup server's root user.

        
ssh xli@128.164.54.240 "mkdir -p /home/back_up/Weekly_Backup/${TODAY}"

rsync --compress --archive --verbose --exclude={"dev","proc","sys","tmp","run","mnt","media","lost+found","cloud_research"} \
	--hard-links --human-readable --inplace --numeric-ids --delete \
	--link-dest=${Destination_Path}/${LASTWEEK} \
	-e ssh ${Folder_Backup} ${Destination} \
	> ${LOG}_$(date +%m%d%Y).log

# Un-hash this if you want to remove old backups (older than 182 days)
OLDBACKUP=`date -d "182 days ago" +"%Y%m%d"`   # old backup has to be multiple of backup inteval. (in here is 7 days.)

ssh xli@128.164.54.240 "rm -R /home/back_up/Weekly_Backup/${OLDBACKUP}"
 

# Writes a log of successful updates
echo "BACKUP success-$TODAY" >> ${LOG}_$(date +%m%d%Y).log

# Clean exit
exit 0
