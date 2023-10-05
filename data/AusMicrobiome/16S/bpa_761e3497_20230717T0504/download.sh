#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=day

cd /data/group/frankslab/project/LFlorence/MycorrhizaAusForests/data/AusMicrobiome/16S/bpa_761e3497_20230717T0504

module load cURL/7.83.0-GCCcore-11.3.0

export CKAN_API_TOKEN=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJiODhUX1hUbUxQUmhpd2NwTVp5LWJpd2c3ZmxFeXlnNm81TGZPcE9CSFZCWmIycUV4aVBnZFJ6VmxNaHdrNFprOUxCR01iTnVYcWVTS1dpTSIsImlhdCI6MTY4OTIxNjA2MH0.W8Y3I0HaNKJIv32682bEPuu5HheP4jWf0jdjytrteSk

# download.sh
# Bulk download tool for the Bioplatforms Australia Data Portal
#
# This UNIX shell script was automatically generated.
#

if [ x"$CKAN_API_TOKEN" = "x" ]; then
  if [ x"$CKAN_API_KEY" = "x" ]; then
    echo "Please set the CKAN_API_TOKEN environment variable."
    echo
    echo "You can create your API Token by browsing to:"
    echo "https://data.bioplatforms.com/user/lukelikesdirt"
    echo
    echo "Go to the API Tokens tab, and generate your token."
    echo
    echo "The API token is a long string of letters and digits"
    echo
    echo "To set the environment variable in Linux/MacOS/Unix, use"
    echo "the following command before running download.sh"
    echo "substituting your API token as required:"
    echo
    echo "export CKAN_API_TOKEN=***********************************"
    echo
    echo "You can check if it has been set correctly with the command:"
    echo
    echo "printenv CKAN_API_TOKEN"
    if [ -t 0 ] ; then
       read -p "Press key to continue... (script will exit) " -n1 -s
    fi
    exit 1
  else
    echo "The API key of the format xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx"
    echo "is now obsolete, and should be replaced wih a freshly generated API Token. "
    echo "You can create your API Token by browsing to:"
    echo "https://data.bioplatforms.com/user/lukelikesdirt"
    echo
    echo "Go to the API Tokens tab, and generate your token."
    echo
    echo "The API token is a long string of letters and digits"
    echo
    echo "To set the environment variable in Linux/MacOS/Unix, use"
    echo "the following command before running download.sh"
    echo "substituting your API token as required:"
    echo
    echo "export CKAN_API_TOKEN=***********************************"
    echo
    echo "You can check if it has been set correctly with the command:"
    echo
    echo "printenv CKAN_API_TOKEN"
     if [ -t 0 ] ; then
       read -p "Press key to continue... script will run using the KEY provided " -n1 -s
     fi
  fi
fi



# Check we are being run from a suitable location

if [ ! -f tmp/bpa_761e3497_20230717T0504_urls.txt ]; then
  echo "tmp/bpa_761e3497_20230717T0504_urls.txt not found"
  echo
  echo "Please change to the directory containing the download.sh script"
  exit 1
fi

# Check for required programs

if ! which curl >/dev/null 2>&1; then
  echo "`curl` is not installed. Please install it."
  echo
  echo "On MacOS, it can be installed via HomeBrew (https://brew.sh/)"
  echo "using the command `brew install curl`"
  exit 1
fi

if ! which md5sum >/dev/null 2>&1; then
  echo "`md5sum` is not installed. Please install it."
  echo
  echo "On MacOS, it can be installed via HomeBrew (https://brew.sh/)"
  echo "using the command `brew install md5sha1sum`"
  exit 1
fi

BPA_AGENT="data.bioplatforms.com download.sh/1.0"

CURL=`which curl`

# if on MacOS, favour homebrew curl over system curl
case "$OSTYPE" in
  darwin*)
    HBCURL="/usr/local/opt/curl/bin/curl"
    if [ -f $HBCURL -a -x $HBCURL ] ; then
        echo "Using curl installed via homebrew"
        CURL="$HBCURL"
    fi
    ;;
  *)
    ;;
esac

# Check program versions

# 7.58 required for correct Authorization header support
CURL_VERSION_REQUIRED="7.58"
CURL_VERSION=$($CURL --version | head -1 | awk '{print $2}')

function max()
{
  local m="$1"
  for n in "$@"
  do
    [ "$n" -gt "$m" ] && m="$n"
  done
  echo "$m"
}

# from https://apple.stackexchange.com/a/261863
function compare_versions()
{
  local v1=( $(echo "$1" | tr '.' ' ') )
  local v2=( $(echo "$2" | tr '.' ' ') )
  
  local len="$(max "${#v1[*]}" "${#v2[*]}")"
  for ((i=0; i<len; i++))
  do
    [ "${v1[i]:-0}" -gt "${v2[i]:-0}" ] && return 1
    [ "${v1[i]:-0}" -lt "${v2[i]:-0}" ] && return 2
  done
  return 0
}

compare_versions $CURL_VERSION $CURL_VERSION_REQUIRED
if [ $? -eq 2 ]; then
  echo "Your 'curl' version is outdated."
  echo
  echo "Path was                   : $CURL"
  echo
  echo "Minimum version required is: $CURL_VERSION_REQUIRED"
  echo "Version available is       : $CURL_VERSION"
  exit 1
fi

# Undertake download

echo "Downloading data"
while read URL; do
  echo "Downloading: $URL"
  if [ x"$CKAN_API_TOKEN" != "x" ]; then
      $CURL -O -L -C - -A "$BPA_AGENT" -H "Authorization: $CKAN_API_TOKEN" "$URL"
  elif [ x"$CKAN_API_KEY" != "x" ]; then
      $CURL -O -L -C - -A "$BPA_AGENT" -H "Authorization: $CKAN_API_KEY" "$URL"
  fi  
  if [ $? -ne 0 ] ; then
     echo "Error downloading: $URL"
  fi
done < tmp/bpa_761e3497_20230717T0504_urls.txt

echo "Data download complete. Verifying checksums:"
md5sum -c tmp/bpa_761e3497_20230717T0504_md5sum.txt 2>&1 | tee tmp/md5sum.log