# Installing PharmGOAT
## From conda
The easiest way to install PharmGOAT is through [conda](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html):

```bash
# create a brand new conda environment and install latest PharmGOAT
conda create -n pharmgoat -c bioconda pharmgoat
# OR install latest into current conda environment
conda install pharmgoat
# OR install a specific version into current conda environment
conda install pharmgoat=0.1.0
```

## From GitHub
Conda updates usually lag the GitHub release by a couple days.
Use the following instructions to get the most recent version directly from GitHub:

1. Navigate to the [latest release](https://github.com/PacificBiosciences/PharmGOAT/releases/latest) and download the tarball file (e.g. `pharmgoat-{version}-x86_64-unknown-linux-gnu.tar.gz`).
2. Decompress the tar file.
3. (Optional) Verify the md5 checksum.
4. Test the binary file by running it with the help option (`-h`).
5. Visit the [User guide](./user_guide.md) for details on running PharmGOAT.

### Example with v0.1.0
```bash
# modify this to update the version
VERSION="v0.1.0"
# get the release file
wget https://github.com/PacificBiosciences/PharmGOAT/releases/download/${VERSION}/pharmgoat-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
# decompress the file into folder ${VERSION}
tar -xzvf pharmgoat-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
cd pharmgoat-${VERSION}-x86_64-unknown-linux-gnu
# optional, check the md5 sum
md5sum -c pharmgoat.md5
# execute help instructions
./pharmgoat -h
```