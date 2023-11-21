# MARS Streamlit Application

## Overview
This Streamlit application integrates the functionalities of MARS (Microbial Abundances Retrieved from Sequencing data) and ANT (Automated NCBI Taxonomy) to provide a comprehensive tool for microbiome researchers. It automates the processing of metagenomic datasets, taxonomic name resolution, and mapping to resources like AGORA2, thereby facilitating a streamlined workflow for microbiome data analysis.

## Installation
To run the application on your local machine, you will need to have Python installed, along with the necessary libraries. Follow these steps to get started:

1. Clone the repository: `git clone https://www.github.com/ines.thiele/mars-pipeline`

2. Navigate to the cloned repository's directory: `cd mars-pipeline`

3. Install the required Python packages: `pip install -r requirements.txt`

## Usage
After installation, you can run the application using Streamlit:

The application will open in your default web browser, or you can access it at `http://localhost:8501`.

# Running the Streamlit Application with Docker

## Prerequisites

Before proceeding, ensure that you have Docker installed on your system. You can download and install Docker from Docker's official website.

## Steps to Run the Application

1. Build the Docker Image
First, you need to build a Docker image from your application. Open a terminal or command prompt and navigate to the directory containing your application's Dockerfile.

Run the following command to build the Docker image:

`docker build -t my_streamlit_app .`

This command builds a new Docker image with the name my_streamlit_app based on the instructions in your Dockerfile. It may take some time as it installs all necessary dependencies.

2. Run the Docker Container
After the image is successfully built, run the container using:

`docker run -p 8501:8501 my_streamlit_app`

This command runs the Docker container and maps port 8501 of the container to port 8501 on your host machine.

3. Accessing the Application
Once the container is running, open your web browser and navigate to:

`http://localhost:8501`

You should now see your Streamlit application running.

## Stopping the Container

To stop the running Docker container, you can press `CTRL+C` in the terminal where the container is running.

Alternatively, you can stop the container using the following steps:

Open a new terminal or command prompt.
Run `docker ps` to list all running containers.
Find the `container ID` of your Streamlit app.
Run docker stop <container_id> to stop the container.
Troubleshooting

If the application is not loading, check if Docker is running correctly and if the correct port is exposed and mapped.
If you encounter errors during the build process, verify that your `Dockerfile` and `requirements.txt` are set up correctly.

## Input Data
The application accepts input data in specific formats. For MARS, a feature table and taxonomy table, or an already combined table are required. Please refer to the sample input files provided in the tests/test_files directory:
`feature_table.txt`, `taxonomy.tsv` and `merged.csv`

## Contributing
Contributions to improve the application are welcome. 

## Support
For support, please open an issue in the GitHub repository or contact ines.thiele@universityofgalway.ie.

## Authors
- Tim Hulshof - Initial work on MARS and ANT
- Bram Nap - Contributions 
- Filippo Martinelli - Contributions

## How to Cite
If you use the MARS Streamlit application in your research, please cite it as follows:
