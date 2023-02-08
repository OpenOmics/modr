## Steps for Building Docker Images

Building a seperate Docker image specific for pycoQC because its python dependencies clash with other tools 
required dependencies.  

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f Dockerfile --tag=pycoqc:v0.1.0 .

# Testing, take a peek inside
docker run -ti pycoqc:v0.1.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag pycoqc:v0.1.0 skchronicles/pycoqc:v0.1.0
docker tag pycoqc:v0.1.0 skchronicles/pycoqc         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/pycoqc:v0.1.0
docker push skchronicles/pycoqc:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan pycoqc:v0.1.0
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
