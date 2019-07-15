#!/usr/bin/env bash
cat << "EOF"
__        __   _                            _
\ \      / /__| | ___ ___  _ __ ___   ___  | |_ ___
 \ \ /\ / / _ \ |/ __/ _ \| '_ ` _ \ / _ \ | __/ _ \
  \ V  V /  __/ | (_| (_) | | | | | |  __/ | || (_) |
   \_/\_/ \___|_|\___\___/|_| |_| |_|\___|  \__\___/

 ____        _____ ___  ____   ____
|  _ \ _   _|  ___/ _ \/ ___| / ___|
| |_) | | | | |_ | | | \___ \| |
|  __/| |_| |  _|| |_| |___) | |___
|_|    \__, |_|   \___/|____/ \____|
       |___/
EOF
CWD=$(pwd -P)
SRCDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
function render_template() {
  eval "echo \"$(cat $1)\""
}
printf "Make sure you have already prepared the data and initiated\niraf before doing the following.\n"
printf "The current working directory is $CWD.\n"
printf "The source directory is $SRCDIR.\n"
printf "============================================================\n"
title="Please select the instrument with which you took the spectra."
prompt="Enter number after the prompt:"
options=("BFOSC, Xinglong 2m16 telescope" "YFOSC, Lijiang 2m4 telescope")

echo "$title"
PS3="$prompt "
select opt in "${options[@]}" "Quit"; do

    case "$REPLY" in

    1 ) printf "============================================================\n"
    echo "You picked $opt which is option $REPLY."
    TELESCOPE="XLT"
    render_template $SRCDIR/config/instrconfig.json.temp > ./myfosc.json
    cp $SRCDIR/iraf_scripts/genlistXLT.cl $CWD
    while true; do
    read -p "Do you want to copy the raw files to ./raw?[y/n]" yn
        case $yn in
          [Yy]* ) bash cpxlraw.sh; echo "Raw data is backuped and copied to ./raw."; exit;;
          [Nn]* ) echo "Otherwise you have to do it yourself."; exit;;
          * ) echo "Please answer y or n.";;
        esac
    done
    break;;
    2 ) printf "============================================================\n"
    TELESCOPE="LJT"
    render_template $SRCDIR/config/instrconfig.json.temp > ./myfosc.json
    cp $SRCDIR/iraf_scripts/genlistLJT.cl $CWD
    echo "You picked $opt which is option $REPLY."
    while true; do
    read -p "Do you want to copy the raw files to ./raw?[y/n]" yn
        case $yn in
          [Yy]* ) bash cpfstext.sh; echo "Raw data is backuped and copied to ./raw."; exit;;
          [Nn]* ) echo "Otherwise you have to do it yourself."; exit;;
          * ) echo "Please answer y or n.";;
        esac
    done
    break;;
    3 ) printf "============================================================\n"
    echo "You picked $opt which is option $REPLY.";;

    $(( ${#options[@]}+1 )) ) echo "Goodbye!"; break;;
    *) echo "Invalid option. Try another one.";continue;;

    esac

done

# while opt=$(zenity --title="$title" --text="$prompt" --list \
#                    --column="Options" "${options[@]}"); do
#
#     case "$opt" in
#     "${options[0]}" ) zenity --info --text="You picked $opt, option 1";;
#     "${options[1]}" ) zenity --info --text="You picked $opt, option 2";;
#     "${options[2]}" ) zenity --info --text="You picked $opt, option 3";;
#     *) zenity --error --text="Invalid option. Try another one.";;
#     esac
#
# done

# STATUS_URI="foobar"
# MONITOR_IP="192.168.1.1"
# # This could also be read in via bash arguments.
# # Google "bash getopts" for more information
# # render a template configuration file
# # expand variables + preserve formatting
# # user="Venkatt"
# # referenced inside the template.txt as $user
# # render_template /path/to/template.txt > path/to/configuration_file
# function render_template() {
#   eval "echo \"$(cat $1)\""
# }
# function generate_httpd_conf {
#   echo "#### Creating /tmp/httpd.conf from template ./httpd.conf.tmpl"
#   render_template httpd.conf.tmpl > /tmp/httpd.conf
# }
