![build workflow](https://github.com/sib-swiss/single-cell-training/actions/workflows/docker-image.yml/badge.svg)
![GitHub Release Date](https://img.shields.io/github/release-date/sib-swiss/single-cell-training)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5703106.svg)](https://doi.org/10.5281/zenodo.5703106)
[![License: CC BY-SA 4.0](https://img.shields.io/badge/License-CC_BY--SA_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by-sa/4.0/)

<!-- Matomo Image Tracker-->
<img referrerpolicy="no-referrer-when-downgrade" src="https://matomo.sib.swiss/matomo.php?idsite=220&amp;rec=1" style="border:0" alt="" />
<!-- End Matomo -->

This website is hosted at: https://sib-swiss.github.io/single-cell-training/

Please refer to [issues](https://github.com/sib-swiss/single-cell-training/issues) for improvements/bugs for course material or the website. 

Any contribution to this course material is highly appreciated :+1:. Please have a look at the [CONTRIBUTING.md](CONTRIBUTING.md) file to learn more on how to contribute. 

# Course website single cell transcriptomics

## Authors

- Tania Wyss [ORCiD](https://orcid.org/0000-0003-2641-0895)
- Rachel Marcone-Jeitziner [ORCiD](https://orcid.org/0000-0002-5711-8435)
- Geert van Geest [ORCiD](https://orcid.org/0000-0002-1561-078X)
- Patricia Palagi [ORCiD](https://orcid.org/0000-0001-9062-6303)

## How reuse this material
This website is generated with quarto. To re-build the website, download and install Rstudio and Quarto CLI. Also make sure you have installed the required packages. After that, clone this repository:

```sh
git clone https://github.com/sib-swiss/single-cell-r-training.git
```
### Using renv
This project uses renv to manage R dependencies. To restore the project's environment, run the following commands in R:

```R
install.packages("renv")
renv::restore()
```
### Render website
Open the project in Rstudio, and run in the terminal to render the full website:

```sh
quarto render
```
Or to preview the website:

```sh
quarto preview
```
To publish the website to GitHub Pages, run:

```sh
quarto publish gh-pages
```
### Run with Docker
A Docker image with all the required software is available on Docker Hub.

To run the Docker container, you can use the following command:

```sh
docker run \
--rm \
-p 8787:8787 \
-v $PWD:/home/rstudio \
sibswiss/training-singlecell-rstudio:latest
```
You can also use the script `Docker/run_locally.sh` to run the container.

### Update Docker image
To add or update R packages in the Docker image, you need to:

1. Add the new packages to the `Docker/install_packages.R` script.
2. Copy the `renv.lock` file from the root directory to the `Docker` directory.
3. Rebuild the Docker image.

#### Using GitHub Actions
You can manually trigger a build of the Docker image by going to the Actions tab of this repository, selecting the "Manual build and push" workflow, and clicking "Run workflow".

#### Locally
Alternatively, you can rebuild the image locally:

```sh
cp renv.lock Docker/
cd Docker
docker build -t sibswiss/training-singlecell-rstudio:latest .
```
