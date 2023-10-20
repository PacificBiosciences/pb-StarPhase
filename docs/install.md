# Installing pb-StarPhase
<!-- ## From conda
The easiest way to install pb-StarPhase is through [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html):

```bash
# create a brand new conda environment and install latest pb-StarPhase
conda create -n pbstarphase -c bioconda pbstarphase
# OR install latest into current conda environment
conda install pbstarphase
# OR install a specific version into current conda environment
conda install pbstarphase=0.7.0
``` -->

## From GitHub
Use the following instructions to get the most recent version of pb-StarPhase directly from GitHub:

1. Navigate to the [latest release](https://github.com/PacificBiosciences/pb-StarPhase/releases/latest) and download the tarball file (e.g. `pbstarphase-{version}-x86_64-unknown-linux-gnu.tar.gz`).
2. Decompress the tar file.
3. (Optional) Verify the md5 checksum.
4. Test the binary file by running it with the help option (`-h`).
5. See the [User guide](./user_guide.md) for details on running pb-StarPhase.

### Example with v0.7.2
```bash
# modify this to update the version
VERSION="v0.7.2"
# get the release file
wget https://github.com/PacificBiosciences/pb-StarPhase/releases/download/${VERSION}/pbstarphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
# decompress the file into folder ${VERSION}
tar -xzvf pbstarphase-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
cd pbstarphase-${VERSION}-x86_64-unknown-linux-gnu
# optional, check the md5 sum
md5sum -c pbstarphase.md5
# execute help instructions
./pbstarphase -h
```