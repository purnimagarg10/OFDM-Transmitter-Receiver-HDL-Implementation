# OFDM-Transmitter-Receiver-HDL-Implementation: MATLAB Simulation to HDL and Vivado Deployment

## Project Overview
This project demonstrates the complete workflow of designing and implementing a wireless communication transceiver, starting from high-level system modeling to hardware implementation on FPGA.

The system is first modeled and verified in **MATLAB**. The validated design is then converted into **HDL (Verilog)** using **HDL Workflow Advisor**. Finally, the generated HDL design is synthesized and implemented using **Xilinx Vivado**.

This project demonstrates the **MATLAB → HDL → FPGA design workflow**, which is widely used in modern wireless communication hardware design and FPGA-based signal processing systems.

---

## Design Workflow

MATLAB Model → HDL Workflow Advisor → Generated Verilog → Vivado FPGA Implementation

---

## System Architecture

The system implements a digital communication transceiver consisting of transmitter and receiver chains.

### Transmitter
The transmitter performs the following operations:

1. Input bit generation  
2. Modulation / symbol mapping  
3. Signal processing operations  
4. Parallel-to-serial conversion  
5. Signal transmission  

### Receiver
The receiver performs the inverse operations:

1. Signal reception  
2. Serial-to-parallel conversion  
3. Demodulation  
4. Signal processing  
5. Bit recovery  

---

## HDL Code Generation Using HDL Workflow Advisor

The MATLAB/Simulink design is converted into synthesizable HDL code using **HDL Workflow Advisor**.

HDL Workflow Advisor performs the following steps:

1. Model compatibility checks for HDL generation  
2. Fixed-point conversion and optimization  
3. HDL code generation  
4. Generation of test benches for verification  
5. Integration support for FPGA tools such as Vivado  

The generated HDL files are then used for FPGA synthesis.

---

## Repository Structure


---

## MATLAB Implementation

MATLAB is used to model and simulate the communication system at an algorithmic level.

Key objectives of the MATLAB implementation include:

- Modeling transmitter and receiver blocks  
- Verifying signal processing operations  
- Validating communication system functionality  
- Generating simulation results before hardware implementation  


---

## Vivado Implementation

The generated HDL design is implemented using **Xilinx Vivado**.

Implementation steps:

1. Create a new Vivado project  
2. Import generated HDL files  
3. Apply timing and pin constraints using `.xdc` files  
4. Run synthesis and implementation  
5. Analyze timing reports and resource utilization  
6. Verify simulation waveforms  

---

## Results

The project produces the following outputs:

- MATLAB simulation results verifying system functionality  
- HDL simulation waveforms validating hardware logic  
- Vivado synthesis reports showing FPGA resource usage  
- Timing analysis ensuring correct hardware performance  

Example outputs include:

- Simulation plots  
- HDL waveform results  
- FPGA synthesis and timing reports  

All relevant results are stored in the **images** folder.

---

## Tools and Technologies Used

- MATLAB  
- HDL Workflow Advisor  
- Verilog HDL  
- Xilinx Vivado  
- Digital Signal Processing concepts  
- FPGA design methodology  

---

## Applications

The concepts demonstrated in this project are used in:

- Wireless communication systems  
- FPGA-based digital signal processing  
- Software defined radio (SDR)  
- Hardware acceleration of signal processing algorithms  

---

## How to Run the Project

### MATLAB Simulation

1. Open MATLAB  
2. Navigate to the `matlab_code` folder  
3. Run the simulation script


---

### Vivado Implementation

1. Open **Xilinx Vivado**  
2. Create a new project  
3. Import the HDL files from the `hdl_code` folder  
4. Add FPGA constraints from `constraints.xdc`  
5. Run synthesis and implementation  

---

## Future Improvements

Potential extensions of this project include:

- Implementing higher-order modulation schemes  
- Adding realistic channel models  
- Optimizing FPGA resource utilization  
- Extending the design to MIMO communication systems  
- Real-time testing on FPGA hardware boards  

---

## Author

Purnima Garg

---

## License

This project is open source and available under the MIT License.
