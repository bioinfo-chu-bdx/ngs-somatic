#!/bin/bash

CURRENT_WOR_DIR=$(pwd)
echo "Current Working Directory $CURRENT_WOR_DIR"

AGENT_HOME=$(cd $(dirname $0) && pwd)

cd $AGENT_HOME

JAVA_HOME=${JAVA_HOME:=/apps/jdk}

#
#return 0 ; ver1 < ver2
#return 1 ; ver1 = ver2
#return 2 ; ver1 > ver2
#
vercomp () {
    if [ $1 == $2 ]
    then
        return 1
    fi
    local IFS=.
    local i ver1=($1) ver2=($2)
    # fill empty fields in ver1 with zeros
    for ((i=${#ver1[@]}; i<${#ver2[@]}; i++))
    do
        ver1[i]=0
    done
    for ((i=0; i<${#ver1[@]}; i++))
    do
        if [ -z ${ver2[i]} ]; then
            # fill empty fields in ver2 with zeros
            ver2[i]=0
        fi
        if [ ${ver1[i]} -gt ${ver2[i]} ]; then
            return 2
        fi
        if [ ${ver1[i]} -lt ${ver2[i]} ]; then
            return 0
        fi
    done
    return 1
}

# Define where is the java executable is
if [ -n "$JAVA_HOME" ] && [ -x "$JAVA_HOME/bin/java" ];  then
    #echo "Found java executable in JAVA_HOME: $JAVA_HOME"
    _java="$JAVA_HOME/bin/java"
elif type -p java; then
    #echo "Found java executable in PATH"
    _java=java
else
    echo "No JAVA_HOME defined and no Java found in PATH. Please intstall Java 8 or higher to use AGeNT"
    exit
fi

if [ "$_java" ]; then
    version=$("$_java" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    vercomp "$version" "1.8"
    if [ $? -lt 1 ]; then
        echo version "$version"
        echo  "Java version is less than 1.8. Version 1.8 or higher is needed for AGeNT"
        exit
    fi
fi


"$_java" -jar lib/Agent*.jar "$CURRENT_WOR_DIR" -home:$AGENT_HOME $@
