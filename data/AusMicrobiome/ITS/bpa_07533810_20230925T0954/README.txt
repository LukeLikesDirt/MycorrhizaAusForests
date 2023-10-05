CKAN Bulk Download
------------------

Search of organization: australian-microbiome bpa_07533810_20230925T0954

Bulk download package generated: 2023-09-25T09:54:54.137382+0000
Number of Organizations        : 1
Number of Packages             : 38
Number of Resources            : 567
Total Space required           : 8.89 GiB
Total Size (bytes)             : 9545188820

This archive contains the following files:

download.py:
Python 3 script, which when executed will download the files,
and then checksum them.  This script is cross platform and is supported
on Linux / MacOS / Windows hosts. Requires the `requests` module.

download.ps1:
Windows PowerShell script, which when executed will download the files,
and then checksum them. There are no dependencies other than PowerShell.

download.sh:
UNIX shell script, which when executed will download the files,
and then checksum them. This is supported on any Linux or MacOS/BSD
system, so long as `curl` is installed.

Before running either of these scripts, please set the CKAN_API_TOKEN
environment variable.

You need to create your API TOKEN within the data portal and copy/save it.
You then use this token whenever you wish to download
data from the data portal.

You can create your API Token by browsing to:
https://data.bioplatforms.com/user/lukelikesdirt, and clicking the API Tokens tab.
Enter a name for your token, then click "Create API Token".

The API key of the format:
xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
is now obsolete, and you should use the instructions above to 
create your new API token.
 
To set the environment variable in Linux/MacOS/Unix, use:
export CKAN_API_TOKEN=xxxxxxxxxxxxxxxxxxx

On Microsoft Windows, within Powershell, use:
$env:CKAN_API_TOKEN="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"


organization_metadata folder:
Contains metadata spreadsheets as CSV for all organizations owning the
selected data resources (files).

package_metadata folder:
Contains metadata spreadsheets as CSV for all selected data packages, grouped
by the type of package (schema). Each data package will contain one or more
resources. This metadata is an amalgamation of all metadata, including
sample contextual metadata and processing metadata.

resource_metadata folder:
Contains metadata spreadsheets as CSV for all selected data resources (files).

QUERY.txt:
Text file which contains metadata about the download results and the original
query

tmp folder:
This folder contains files required by the download scripts. Its
contents can be ignored.


Note all CSV files are encoded as UTF-8 with a Byte Order Mark (BOM) to
enable character set detection by recent versions of Microsoft Excel.
