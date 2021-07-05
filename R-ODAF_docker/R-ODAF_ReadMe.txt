
# INSTALL:
# In order to use the supplied dockerfile, docker must be installed on your system:
# A) Depending on your system, you need to put sudo before each command. Since its not always nessesary, its written here between ()
(sudo) apt install docker.io
# B) Docker needs to be setuo to run at startup:
(sudo) systemctl start docker
(sudo) systemctl enable docker
# C.1) Thereafter the docker environment needs to be build:
(sudo) docker build -t rodaf_sequencing /Path_to_DOCKERFILE/
# C.2) or, if run from the location where the dockerfile is located:
(sudo) docker build -t rodaf_sequencing .

# RUNNING THE DOCKER:
# Step 1: Start the docker
(sudo) docker run -v /Put_Here_The_Work_Directory_On_Your_System:/data -it rodaf_sequencing

# Step 2: Make sure all that the next files are located within the folder to which the docker is mounted (/Put_Here_The_Work_Directory_On_Your_System):
         # sample fastq files
         # genomefiles (fasta + gtf)
         # script used by the docker (R-ODAF_sequencing_DataPreprocessing.sh)

# Step 3: Adapt the top section of the R-ODAF_sequencing_DataPreprocessing.sh (parameters & filelocations)

# Step 4: Make sure all files and folders are accessible (if not, run the next line in the /Put_Here_The_Work_Directory_On_Your_System folder:
(sudo) chmod -R 777 *

# Step 5: Run the R-ODAF_sequencing_DataPreprocessing.sh
cd data
./R-ODAF_sequencing_DataPreprocessing.sh
# If the docker stops due to permission errors, repeat step 4

# Step 6: After the docker has finished, the permissions need to be reset in order to access the output:
(sudo) chmod -R 777 *

# Step 7: proceed with Differential expression analysis. An example script for this purpose has been supplied on Github: R-ODAFsequencing_DEGidentification.R 

### Other usefull commands
# Check created dockers:
sudo docker container ls -a
# Remove a finisched docker:
sudo docker container rm Put_Docker_ID_here






