## HCA

Separate project folder for the HCA (Human Chromosome Analysis) and Bacterial chromosome analysis tool


### How to Run the HCA Barcode Alignment Tool with Path Inclusion:

1. **Download the Code**:
   - Ensure you have downloaded the HCA Barcode Alignment tool (https://github.com/dnadevcode/hca) code to your local machine.

2. **Open MATLAB**:
   - Launch MATLAB on your computer.

3. **Navigate to Code Directory**:
   - Use the MATLAB command window or the file explorer to navigate to the directory where you have saved the HCA Barcode Alignment tool code.

4. **Add Paths**:
   - Use the `addpath` function to include the necessary directories. For example:
    addpath(genpath(pwd));


## Tutorial
### How to Use GUI Buttons in HCA Barcode Alignment Tool:

1. **Launch the HCA Barcode Alignment Tool**:
   - Execute the provided MATLAB function (hca_barcode_alignment) to launch the HCA Barcode Alignment tool. You can do this by running the function in MATLAB.

2. **Navigate Through the Menu**:
   - Access the HCA menu by clicking on its corresponding button in the GUI interface.

#### Theory Menu:
   - Under the HCA menu, you'll find options related to theory generation and manipulation.

   a. **Generate Theory**:

      - Click on the "Theory" submenu.
      - Choose "Generate theory" to initiate the theory generation process.
      - Provide necessary inputs such as file paths or settings.
      - Click on the "Run" button to start theory generation.
      - Monitor progress and review the results.

   b. **Concatenate Theories**:

        - **Click the "Run Concatenate Theories" Button**: Locate and click on the "Run Concatenate Theories" button in your MATLAB user interface.

        - **Select Files**: After clicking the button, a dialog box will appear prompting you to select one or more files. These files should contain theory data that you want to merge together.

        - **Wait for Processing**: Once you've selected the files, the function will start processing them. This may take a moment depending on the size and number of files.

        - **Review Messages**: As the function runs, it will display messages in the MATLAB command window. If all selected files have compatible settings, the function will concatenate the theories. If there are any issues, such as incompatible settings between files, the function will inform you about them.

        - **Save Concatenated Data**: If the function successfully merges the theories, it will prompt you to save the concatenated data. Choose a location and filename for the new file.

        - **Done!**: Once the data is saved, the function will confirm the successful concatenation. You can then access the merged theory data from the file you saved it to.


#### Alignment Menu:
   - Explore alignment-related options and functionalities.

   a. **Run Alignment**:

      - Access the "Alignment" submenu.
      - Choose "Run Alignment" to initiate alignment.
      - Provide required inputs and follow on-screen prompts.
      - Monitor progress and review the alignment results.

   b. **Shrink filter**:

      - Under "Alignment", select "Shrink filter" to begin shrinking filtering process.
      - Follow instructions and provide necessary inputs.
      - Monitor progress and review the results.

   c. **Duplicates Sorter**:

      - Choose "Duplicates sorter" to start duplicates sorting.
      - Provide inputs as needed and follow on-screen instructions.
      - Review the sorting results once completed.

   d. **Save Alignment Result**:

      - Select "Save Alignment Result" to save the alignment session.
      - Specify a location and filename to store the session data.

   e. **Load Session**:

      - Use "Load session" to open a saved alignment session.
      - Navigate to the saved session file and load it for further analysis.

### Note:

- Follow the provided instructions within the GUI interface for each button to ensure correct usage.
- Make sure you have the required permissions and access to the necessary data files.
- Consult the tool's documentation or seek assistance for any issues or questions regarding specific functionalities.