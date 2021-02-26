# OASIS Submission - Island Scanning

## Build & Run Instructions
The following instructions should get you up and running:
1. Build and run the Dockerfile via `docker build -t osu_cdme:latest -m 2GB .`
2. Run the Dockerfile via `docker run -it osu_cdme`
3. Get this directory into the container. Options include `docker cp` and `git clone`.
4. Build the source code by running `build.bat`
5. Execute layer/scan generation for the build plate by running `run_Build_Plate.bat` and for the NIST Plate by running `run_NIST.bat`
6. You can then find the output .scn file in the location specified in the respective config file.

## Scanpath Overview
This submission is an implementation of the well-researched island scanpath algorithm. The algorithm splits vectors along evenly spaced horizontal and vertical lines, creating squares, then iterates through the vectors in those squares in order. 