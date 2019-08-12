## Low CBF Model

This repository contains a Matlab model of signal processing performed by Low-CBF FPGA firmware. The model can be run using either Matlab or Octave.

## About

Low-CBF FPGAs receive UDP data packets containing station data from LFAA stations. They filter, beamform, and cross-correlate the data and send output to SDP. The model:
* Simulates the sky and LFAA stations to generate synthetic data packets that LFAA stations would produce
* Processes the packets exactly as the FPGA implemenation of Low-CBF does
* Generates outputs and intermediate results that can be compared with real FPGA outputs

There is a related tool (lfaa-sim) that is able to take packet output from this model and play the packets over a 40GbE network. This allows an FPGA implementing Low-CBF processing to be tested in real-time by capturing the outputs that occur in response to specific simulated circumstances and comparing them with the model's intermediate results and outputs. See https://github.com/ska-telescope/low-cbf-tools

## Documentation
Documentation of the model, and a 'howto' describing its use can be found here: https://confluence.skatelescope.org/display/SE/Packetised+Model+Overview
