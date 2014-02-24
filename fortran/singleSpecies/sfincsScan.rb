#!/usr/bin/env ruby

$filename="input.namelist"
jobFilename="job.sfincsScan"

require 'fileutils'
include FileUtils        
            
if !File.exists?($filename)
   puts "Error! #{$filename} not found."
   exit
end
puts "File #{$filename} exists."

if !File.exists?(jobFilename)
   puts "Error! #{$jobFilename} not found."
   exit
end
puts "File #{jobFilename} exists."

# Make sure job file does not already contain "#PBS -N"
inFile = File.open(jobFilename,"r")
lines = inFile.readlines
lines.each {|line| 
   if (line.include? "#PBS") & (line.include? "-N")
      puts "Error! #{jobFilename} should not include a -N line with job name."
      exit
   end
}
inFile.close

def namelistLineContains(line, variableName)
   line = line.strip.downcase
   variableName = variableName.downcase

   if (!line.include? variableName) : return false end
   if (line[0].chr == "!") : return false end

   nextChar = line[variableName.size].chr
   if (line[0..(variableName.size-1)] == variableName) & ((nextChar==" ") | (nextChar == "="))
      return true
   else
      return false
   end
end

def readInput(variableName, intOrFloat)
   # set intOrFloat = 0 to read an integer or 1 to read a float.

   if !variableName.is_a?(String)
      puts "Error: variableName must be a string."
      exit
   end

   if !intOrFloat.is_a?(Integer)
      puts "Error: intOrFloat must be an integer."
      exit
   end
   if (intOrFloat > 1) | (intOrFloat < 0)
      puts "Error: intOrFloat must be 0 or 1."
      exit
   end

   variableName = variableName.downcase
   s = `grep -i #{variableName} #{$filename}`
   numMatches = 0
   value = 0
   s.each_line do |line|
      line = line.strip.downcase
      # If not a comment:
      if line[0].chr != "!"
         nextChar = line[variableName.size].chr
         if (line[0..(variableName.size-1)] == variableName) & ((nextChar==" ") | (nextChar == "="))
            # Ruby doesn't understand the fortran scientific notation 1d-0 so replace it with 1e-0.
            substring = line[(line.index("=")+1)..(line.size-1)]
            substring.gsub!("d-","e-")
            #value = line.scan(/\d+/)[0].to_f
            if intOrFloat==0
               value = substring.to_i
            else
               value = substring.to_f
            end
            numMatches = numMatches + 1
         end
      end
   end
   if numMatches < 1
      puts "Error! No lines in #{$filename} match #{variableName}."
      exit
   end
   if numMatches > 1
      puts "Warning: more than 1 line in #{$filename} matches #{variableName}."
   end
   if intOrFloat == 0
      puts "Read #{variableName} = " + value.to_i.to_s
      return value.to_i
   else
      puts "Read #{variableName} = " + value.to_s
      return value
   end
end

def linspace(min, max, nn)
   return (0..(nn-1)).collect{|x| x*(max-min)/(nn-1.0)+min}
end

def logspace(min, max, nn)
   if (min <= 0)
      puts "Error! In logspace, min must be positive"
      exit
   end
   if (max <= 0)
      puts "Error! In logspace, max must be positive"
      exit
   end
   logs = linspace(Math.log(min), Math.log(max), nn)
   return logs.collect {|x| Math.exp(x)}
end

programMode = readInput("programMode",0)

case programMode
when 2
   puts "Beginning convergence scan"

   Ntheta = readInput("Ntheta",0)
   NthetaMinFactor = readInput("NthetaMinFactor",1)
   NthetaMaxFactor = readInput("NthetaMaxFactor",1)
   NthetaNumRuns = readInput("NthetaNumRuns",0)
   Nthetas_tmp = logspace(Ntheta*NthetaMinFactor, Ntheta*NthetaMaxFactor, NthetaNumRuns).collect{|x| x.round}
   # Force Ntheta to be odd:
   Nthetas = Nthetas_tmp.collect{|x| if (x % 2 == 1) then x else x+1 end}
   puts "Nthetas:"
   p Nthetas

   Nzeta = readInput("Nzeta",0)
   NzetaMinFactor = readInput("NzetaMinFactor",1)
   NzetaMaxFactor = readInput("NzetaMaxFactor",1)
   NzetaNumRuns = readInput("NzetaNumRuns",0)
   Nzetas_tmp = logspace(Nzeta*NzetaMinFactor, Nzeta*NzetaMaxFactor, NzetaNumRuns).collect{|x| x.round}
   # Force Nzeta to be odd:
   Nzetas = Nzetas_tmp.collect{|x| if (x % 2 == 1) then x else x+1 end}
   puts "Nzetas:"
   p Nzetas

   Nxi = readInput("Nxi",0)
   NxiMinFactor = readInput("NxiMinFactor",1)
   NxiMaxFactor = readInput("NxiMaxFactor",1)
   NxiNumRuns = readInput("NxiNumRuns",0)
   Nxis = logspace(Nxi*NxiMinFactor, Nxi*NxiMaxFactor, NxiNumRuns).collect{|x| x.round}
   puts "Nxis:"
   p Nxis

   Nx = readInput("Nx",0)
   NxMinFactor = readInput("NxMinFactor",1)
   NxMaxFactor = readInput("NxMaxFactor",1)
   NxNumRuns = readInput("NxNumRuns",0)
   Nxs = logspace(Nx*NxMinFactor, Nx*NxMaxFactor, NxNumRuns).collect{|x| x.round}
   puts "Nxs:"
   p Nxs

   NL = readInput("NL",0)
   NLMinFactor = readInput("NLMinFactor",1)
   NLMaxFactor = readInput("NLMaxFactor",1)
   NLNumRuns = readInput("NLNumRuns",0)
   NLs = logspace(NL*NLMinFactor, NL*NLMaxFactor, NLNumRuns).collect{|x| x.round}
   puts "NLs:"
   p NLs

   NxPotentialsPerVth = readInput("NxPotentialsPerVth",1)
   NxPotentialsPerVthMinFactor = readInput("NxPotentialsPerVthMinFactor",1)
   NxPotentialsPerVthMaxFactor = readInput("NxPotentialsPerVthMaxFactor",1)
   NxPotentialsPerVthNumRuns = readInput("NxPotentialsPerVthNumRuns",0)
   NxPotentialsPerVths = logspace(NxPotentialsPerVth*NxPotentialsPerVthMinFactor, NxPotentialsPerVth*NxPotentialsPerVthMaxFactor, NxPotentialsPerVthNumRuns)
   puts "NxPotentialsPerVths:"
   p NxPotentialsPerVths

   xMax = readInput("xMax",1)
   xMaxMinFactor = readInput("xMaxMinFactor",1)
   xMaxMaxFactor = readInput("xMaxMaxFactor",1)
   xMaxNumRuns = readInput("xMaxNumRuns",0)
   xMaxs = logspace(xMax*xMaxMinFactor, xMax*xMaxMaxFactor, xMaxNumRuns)
   puts "xMaxs:"
   p xMaxs

   solverTolerance = readInput("solverTolerance",1)
   solverToleranceMinFactor = readInput("solverToleranceMinFactor",1)
   solverToleranceMaxFactor = readInput("solverToleranceMaxFactor",1)
   solverToleranceNumRuns = readInput("solverToleranceNumRuns",0)
   solverTolerances = logspace(solverTolerance*solverToleranceMinFactor, solverTolerance*solverToleranceMaxFactor, solverToleranceNumRuns)
   puts "solverTolerances:"
   p solverTolerances


   numRunsInScan = 1 + NthetaNumRuns + NzetaNumRuns + NxiNumRuns \
       + NLNumRuns + NxNumRuns + NxPotentialsPerVthNumRuns + xMaxNumRuns + solverToleranceNumRuns

   baseCase = [Ntheta,Nzeta,Nxi,NL,Nx,NxPotentialsPerVth,xMax,solverTolerance];
   parametersForScan = Array.new
   for i in 1..numRunsInScan
      # Need to use .dup since otherwise each element of the outer array will point to the same instance of the inner array.
      parametersForScan << baseCase.dup
   end

   currentIndex = 1

   for i in 1..NthetaNumRuns
      parametersForScan[currentIndex][0] = Nthetas[i-1]
      currentIndex += 1
   end

   for i in 1..NzetaNumRuns
      parametersForScan[currentIndex][1] = Nzetas[i-1]
      currentIndex += 1
   end

   for i in 1..NxiNumRuns
      parametersForScan[currentIndex][2] = Nxis[i-1]
      currentIndex += 1
   end

   for i in 1..NLNumRuns
      parametersForScan[currentIndex][3] = NLs[i-1]
      currentIndex += 1
   end

   for i in 1..NxNumRuns
      parametersForScan[currentIndex][4] = Nxs[i-1]
      currentIndex += 1
   end

   for i in 1..NxPotentialsPerVthNumRuns
      parametersForScan[currentIndex][5] = NxPotentialsPerVths[i-1]
      currentIndex += 1
   end

   for i in 1..xMaxNumRuns
      parametersForScan[currentIndex][6] = xMaxs[i-1]
      currentIndex += 1
   end

   for i in 1..solverToleranceNumRuns
      parametersForScan[currentIndex][7] = solverTolerances[i-1]
      currentIndex += 1
   end

   if (currentIndex != numRunsInScan)
      puts "Error! Something went wrong."
      exit
   end

   # Now remove any duplicates

   i = 0
   while i < numRunsInScan-1
      j=i+1
      while j < numRunsInScan
         if parametersForScan[i] == parametersForScan[j]
            # Item j is a duplicate, so remove it.
            parametersForScan.delete_at(j)
            numRunsInScan -= 1
	    j -= 1
         end
         j += 1
      end
      i += 1
   end

else
   puts "I do not know what to do with programMode = "+programMode.to_s
   exit
end

puts "Parameters for scan:"
p parametersForScan

puts "Scan will consist of #{numRunsInScan} runs."

#exit

for i in 0..(numRunsInScan-1)
   dirName = i.to_s
   mkdir(dirName)
   
   # Copy pbs file
   outFilename = dirName + "/" + jobFilename
   inFile = File.open(jobFilename,"r")
   outFile = File.open(outFilename,"w")
   lines = inFile.readlines
   outFile.write(lines[0])
   outFile.write("#PBS -N sfincs.#{i}\n")
   for j in 1..(lines.size-1)
      outFile.write(lines[j])
   end
   inFile.close
   outFile.close

   # Copy input.namelist
   outFilename = dirName + "/" + $filename
   inFile = File.open($filename,"r")
   outFile = File.open(outFilename,"w")
   lines = inFile.readlines
   for j in 0..(lines.size-1)
      line = lines[j].strip

      if namelistLineContains(line,"programMode")
         line = "programMode = 1"
      end

      if namelistLineContains(line,"Ntheta")
         line = "Ntheta = " + parametersForScan[i][0].to_s
      end

      if namelistLineContains(line,"Nzeta")
         line = "Nzeta = " + parametersForScan[i][1].to_s
      end

      if namelistLineContains(line,"Nxi")
         line = "Nxi = " + parametersForScan[i][2].to_s
      end

      if namelistLineContains(line,"NL")
         line = "NL = " + parametersForScan[i][3].to_s
      end

      if namelistLineContains(line,"Nx")
         line = "Nx = " + parametersForScan[i][4].to_s
      end

      if namelistLineContains(line,"NxPotentialsPerVth")
         line = "NxPotentialsPerVth = " + parametersForScan[i][5].to_s
      end

      if namelistLineContains(line,"xMax")
         line = "xMax = " + parametersForScan[i][6].to_s
      end

      if namelistLineContains(line,"solverTolerance")
         line = "solverTolerance = " + parametersForScan[i][7].to_s
      end

      outFile.write(line + "\n")
   end
   inFile.close
   outFile.close

   # Submit job!
   puts "Submitting job #{i}"
   job1 = fork do
      exec "cd #{i}; qsub #{jobFilename}"
   end
   Process.detach(job1)
   #puts `cd #{i}; qsub #{jobFilename} &`
   puts "Done submitting job #{i}"
end

puts "Finished submitting jobs for scan."