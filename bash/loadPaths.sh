dep_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )/../dep/ ## get current dir, then dependency dir

## Load dependency paths
java=${dep_dir}/jre1.8.0_191/bin/java
snpsift=${dep_dir}/snpEff/SnpSift.jar