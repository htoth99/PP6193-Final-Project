## Submission Notes

Hi Jelmer! This note is to help you understand everything that's in my project directory. I have instructions below and then some notes for you below the instructions. Some parts of the README.md and USAGE.md convey the same information, but will change after the project to include separate information. 

1. Read submission note
2. Read README.md again (changed it a bit from proposal)
3. Read USAGE.md for instructions to run the entire workflow
4. Create envrionments (detailed in USAGE.md) and copy data from location noted below
5. Run workflow!

## Other important notes
For raw data, copy from this location
```bash
cp /fs/ess/PAS2700/users/htoth99/PP6193_FinalProject/rawdata
```
I moved the project structure to the USAGE.md beacuse I felt that it fit there better.

I still haven't decided which binner I like the best. I need to do host removal to really understand which one would be best for plant pathogenic work. Therefore, I've kept both of them in here because I like both. There's no read mapping for maxbin2, which is super nice, but it does take more memory/time in the end. Since I am indecisive, you'll find two directories within binning for metabat2 and maxbin2 to separate the outputs of both. 

### What doesn't this workflow include? 
This workflow also does not include extensive error messaging/guiding. This is something I will continue to work on once I've established and utilized this pipeline across many samples. I also want to have others test it this out so I can see the error others may be encountering. Since the set up currently works for me, and is based on copying and pasting from the runner scripts, it was challenging for me to generate errors! Definitely something I want to introduce once I start handing this off to people.

I'm planning on adding checkm and host removal to this workflow as well, but I need to figure out how to do host removal first. Coming soon for sure!
