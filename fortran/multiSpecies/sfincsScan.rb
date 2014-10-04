#!/usr/bin/env ruby

# This ruby script launches a bunch of sfincs jobs in a batch system,
# with the capability to perform several types of scans.

# Scan types (i.e. values for programMode) implemented in this script:
# 2 = Convergence scan
# 8 = Scan of dPhiHatdpsiN, at fixed resolutions.
# 9 = Simultaneous convergence scan and scan of dPhiHatdpsiN
# 10 = Simultaneous scan of dphiHatdpsi, of collisionOperator=0 and 1, and of DKES vs full trajectories, keeping resolutions fixed.

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
  
  if (!line.include? variableName)
    return false 
  end
  if (line[0].chr == "!") 
    return false 
  end
  
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
    raise "Error! In logspace, min must be positive"
  end
  if (max <= 0)
    raise "Error! In logspace, max must be positive"
  end
  logs = linspace(Math.log(min), Math.log(max), nn)
  return logs.collect {|x| Math.exp(x)}
end

programMode = readInput("programMode",0)

case programMode
when 2,9
  puts "Beginning convergence scan"
  
  Ntheta = readInput("Ntheta",0)
  NthetaMinFactor = readInput("NthetaMinFactor",1)
  NthetaMaxFactor = readInput("NthetaMaxFactor",1)
  NthetaNumRuns0 = readInput("NthetaNumRuns",0)
  Nthetas_tmp = logspace(Ntheta*NthetaMinFactor, Ntheta*NthetaMaxFactor, NthetaNumRuns0).collect{|x| x.round}
  # Force Ntheta to be odd:
  Nthetas = Nthetas_tmp.collect{|x| if (x % 2 == 1) then x else x+1 end}.uniq
  NthetaNumRuns = Nthetas.size
  puts "Nthetas:"
  p Nthetas
  
  Nzeta = readInput("Nzeta",0)
  NzetaMinFactor = readInput("NzetaMinFactor",1)
  NzetaMaxFactor = readInput("NzetaMaxFactor",1)
  NzetaNumRuns0 = readInput("NzetaNumRuns",0)
  Nzetas_tmp = logspace(Nzeta*NzetaMinFactor, Nzeta*NzetaMaxFactor, NzetaNumRuns0).collect{|x| x.round}
  # Force Nzeta to be odd:
  Nzetas = Nzetas_tmp.collect{|x| if (x % 2 == 1) then x else x+1 end}.uniq
  NzetaNumRuns = Nzetas.size
  puts "Nzetas:"
  p Nzetas
  
  Nxi = readInput("Nxi",0)
  NxiMinFactor = readInput("NxiMinFactor",1)
  NxiMaxFactor = readInput("NxiMaxFactor",1)
  NxiNumRuns0 = readInput("NxiNumRuns",0)
  Nxis = logspace(Nxi*NxiMinFactor, Nxi*NxiMaxFactor, NxiNumRuns0).collect{|x| x.round}.uniq
  NxiNumRuns = Nxis.size
  puts "Nxis:"
  p Nxis
  
  Nx = readInput("Nx",0)
  NxMinFactor = readInput("NxMinFactor",1)
  NxMaxFactor = readInput("NxMaxFactor",1)
  NxNumRuns0 = readInput("NxNumRuns",0)
  Nxs = logspace(Nx*NxMinFactor, Nx*NxMaxFactor, NxNumRuns0).collect{|x| x.round}.uniq
  NxNumRuns = Nxs.size
  puts "Nxs:"
  p Nxs
  
  NL = readInput("NL",0)
  NLMinFactor = readInput("NLMinFactor",1)
  NLMaxFactor = readInput("NLMaxFactor",1)
  NLNumRuns0 = readInput("NLNumRuns",0)
  NLs = logspace(NL*NLMinFactor, NL*NLMaxFactor, NLNumRuns0).collect{|x| x.round}.uniq
  NLNumRuns = NLs.size
  puts "NLs:"
  p NLs
  
  NxPotentialsPerVth = readInput("NxPotentialsPerVth",1)
  NxPotentialsPerVthMinFactor = readInput("NxPotentialsPerVthMinFactor",1)
  NxPotentialsPerVthMaxFactor = readInput("NxPotentialsPerVthMaxFactor",1)
  NxPotentialsPerVthNumRuns0 = readInput("NxPotentialsPerVthNumRuns",0)
  NxPotentialsPerVths = logspace(NxPotentialsPerVth*NxPotentialsPerVthMinFactor, NxPotentialsPerVth*NxPotentialsPerVthMaxFactor, NxPotentialsPerVthNumRuns0).uniq
  NxPotentialsPerVthNumRuns = NxPotentialsPerVths.size
  puts "NxPotentialsPerVths:"
  p NxPotentialsPerVths
  
  xMax = readInput("xMax",1)
  xMaxMinFactor = readInput("xMaxMinFactor",1)
  xMaxMaxFactor = readInput("xMaxMaxFactor",1)
  xMaxNumRuns0 = readInput("xMaxNumRuns",0)
  xMaxs = logspace(xMax*xMaxMinFactor, xMax*xMaxMaxFactor, xMaxNumRuns0).uniq
  xMaxNumRuns = xMaxs.size
  puts "xMaxs:"
  p xMaxs
  
  solverTolerance = readInput("solverTolerance",1)
  solverToleranceMinFactor = readInput("solverToleranceMinFactor",1)
  solverToleranceMaxFactor = readInput("solverToleranceMaxFactor",1)
  solverToleranceNumRuns0 = readInput("solverToleranceNumRuns",0)
  solverTolerances = logspace(solverTolerance*solverToleranceMinFactor, solverTolerance*solverToleranceMaxFactor, solverToleranceNumRuns0).uniq
  solverToleranceNumRuns = solverTolerances.size
  puts "solverTolerances:"
  p solverTolerances
  
  
  numRunsInScan = 1 + NthetaNumRuns + NzetaNumRuns + NxiNumRuns \
  + NLNumRuns + NxNumRuns + NxPotentialsPerVthNumRuns + xMaxNumRuns + solverToleranceNumRuns
  
  baseCase = [Ntheta,Nzeta,Nxi,NL,Nx,NxPotentialsPerVth,xMax,solverTolerance];
  parametersForScan = Array.new
  descriptions = Array.new
  for i in 1..numRunsInScan
    # Need to use .dup since otherwise each element of the outer array will point to the same instance of the inner array.
    parametersForScan << baseCase.dup
  end
  
  currentIndex = 1
  descriptions[0]="baseCase"
  
  for i in 1..NthetaNumRuns
    parametersForScan[currentIndex][0] = Nthetas[i-1]
    descriptions[currentIndex] = "Ntheta" + Nthetas[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NzetaNumRuns
    parametersForScan[currentIndex][1] = Nzetas[i-1]
    descriptions[currentIndex] = "Nzeta" + Nzetas[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NxiNumRuns
    parametersForScan[currentIndex][2] = Nxis[i-1]
    descriptions[currentIndex] = "Nxi" + Nxis[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NLNumRuns
    parametersForScan[currentIndex][3] = NLs[i-1]
    descriptions[currentIndex] = "NL" + NLs[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NxNumRuns
    parametersForScan[currentIndex][4] = Nxs[i-1]
    descriptions[currentIndex] = "Nx" + Nxs[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..NxPotentialsPerVthNumRuns
    parametersForScan[currentIndex][5] = NxPotentialsPerVths[i-1]
    descriptions[currentIndex] = "NxPotentialsPerVth" + NxPotentialsPerVths[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..xMaxNumRuns
    parametersForScan[currentIndex][6] = xMaxs[i-1]
    descriptions[currentIndex] = "xMax" + xMaxs[i-1].to_s
    currentIndex += 1
  end
  
  for i in 1..solverToleranceNumRuns
    parametersForScan[currentIndex][7] = solverTolerances[i-1]
    descriptions[currentIndex] = "solverTolerance" + solverTolerances[i-1].to_s
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
      isThisADuplicate = true
      for k in 0...(parametersForScan[0].size)
        # If values are within 0.1%, the difference is probably just due to roundoff error,
        # so treat the values as the same:
        if (parametersForScan[i][k] - parametersForScan[j][k]*1.0).abs/parametersForScan[i][k] > 0.001
          isThisADuplicate = false
        end
      end
      if isThisADuplicate
        # Item j is a duplicate, so remove it.
        parametersForScan.delete_at(j)
        descriptions.delete_at(j)
        numRunsInScan -= 1
        j -= 1
      end
      j += 1
    end
    i += 1
  end
  
  puts "Parameters for scan:"
  p parametersForScan
  
  puts "Description of runs:"
  p descriptions
  
when 8
  # Scan of E_r (i.e. of dPhiHatdpsiN)
  
  NErs = readInput("NErs",0)
  dPhiHatdpsiN_min = readInput("dPhiHatdpsiN_min",1)
  dPhiHatdpsiN_max = readInput("dPhiHatdpsiN_max",1)
  dPhiHatdpsiNs = linspace(dPhiHatdpsiN_min, dPhiHatdpsiN_max, NErs)
  
  numRunsInScan = NErs
  
  puts "dPhiHatdpsiNs:"
  p dPhiHatdpsiNs

when 10
  descriptions = ["FP_full", "FP_DKES", "PAS_full", "PAS_DKES"]
  numRunsInScan = 4
  
else
  puts "I do not know what to do with programMode = "+programMode.to_s
  exit
end


case programMode
when 9,10
  # A 2-level scan, with Er as the inner directory.
  
  NErs = readInput("NErs",0)
  dPhiHatdpsiN_min = readInput("dPhiHatdpsiN_min",1)
  dPhiHatdpsiN_max = readInput("dPhiHatdpsiN_max",1)
  dPhiHatdpsiNs = linspace(dPhiHatdpsiN_min, dPhiHatdpsiN_max, NErs)
  
  puts "dPhiHatdpsiNs:"
  p dPhiHatdpsiNs
  
  puts "Scan will consist of #{numRunsInScan} runs for each of #{NErs} Ers, for a total of #{numRunsInScan*NErs} runs."
  
  for i in 0...numRunsInScan
    outerDirName = descriptions[i]
    if !File.exists?(outerDirName)
      mkdir(outerDirName)
    end
    
    puts "Beginning to submit jobs for "+descriptions[i]
    
    for k in 0...NErs
      innerDirNum = k-1
      # Starting with 0, look for the first unused whole number to use for a directory name:
      begin
        innerDirNum += 1
      	innerDirName = innerDirNum.to_s
        dirName = outerDirName + "/" + innerDirName
      end until !File.exists?(dirName)
      mkdir dirName
      
      # Copy pbs file
      outFilename = dirName + "/" + jobFilename
      inFile = File.open(jobFilename,"r")
      outFile = File.open(outFilename,"w")
      lines = inFile.readlines
      outFile.write(lines[0])
      outFile.write("#PBS -N #{outerDirName}.#{innerDirNum}\n")
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

        case programMode
        when 9
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

        when 10
          if i==0 or i==1
            # Use Fokker-Planck collisions:
            if namelistLineContains(line,"collisionOperator")
              line = "collisionOperator = 0"
            end
          else
            # Use pure pitch-angle scattering collisions:
            if namelistLineContains(line,"collisionOperator")
              line = "collisionOperator = 1"
            end
          end

          if i==0 or i==2
            # Use full trajectories:

            if namelistLineContains(line,"includeXDotTerm")
              line = "includeXDotTerm = .true."
            end
            if namelistLineContains(line,"includeElectricFieldTermInXiDot")
              line = "includeElectricFieldTermInXiDot = .true."
            end
            if namelistLineContains(line,"useDKESExBDrift")
              line = "useDKESExBDrift = .false."
            end
            if namelistLineContains(line,"include_fDivVE_term")
              line = "include_fDivVE_term = .false."
            end

          else
            # Use DKES trajectories:

            if namelistLineContains(line,"includeXDotTerm")
              line = "includeXDotTerm = .false."
            end
            if namelistLineContains(line,"includeElectricFieldTermInXiDot")
              line = "includeElectricFieldTermInXiDot = .false."
            end
            if namelistLineContains(line,"useDKESExBDrift")
              line = "useDKESExBDrift = .true."
            end
            if namelistLineContains(line,"include_fDivVE_term")
              line = "include_fDivVE_term = .false."
            end

          end

        else
          raise "Program should not get here."
        end

        if namelistLineContains(line,"dPhiHatdpsiN")
          line = "dPhiHatdpsiN = " + dPhiHatdpsiNs[k].to_s
        end
        
        
        outFile.write(line + "\n")
      end
      inFile.close
      outFile.close
      
      # Submit job!
      puts "Submitting job #{innerDirNum} in the Er scan."
      # Make sure base case jobs get submitted first:
      #if i==0 and programMode == 9
      if true
        # In this way of submitting a job, ruby waits for qsub to complete before moving on.
        puts `cd #{dirName}; qsub #{jobFilename} &`
      else
        # In this way of submitting a job, ruby does not wait for qsub to complete before moving on.
        job1 = fork do
          exec "cd #{dirName}; qsub #{jobFilename}"
        end
        Process.detach(job1)
      end
      puts "Done submitting job #{innerDirNum}"
    end # of loop j over Ers
  end # of loop i over convergence scan
  
else
  # A standard scan, like programMode = 2 or 8 but not 9 or 10.
  
  puts "Scan will consist of #{numRunsInScan} runs."
  
  
  for i in 0..(numRunsInScan-1)
    dirNum = i-1
    # Starting with 0, look for the first unused whole number to use for a directory name:
    begin
      dirNum += 1
      dirName = dirNum.to_s
    end until !File.exists?(dirName)
    mkdir(dirName)
    
    # Copy pbs file
    outFilename = dirName + "/" + jobFilename
    inFile = File.open(jobFilename,"r")
    outFile = File.open(outFilename,"w")
    lines = inFile.readlines
    outFile.write(lines[0])
    outFile.write("#PBS -N sfincs.#{dirNum}\n")
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
      
      case programMode
      when 2
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
        
      when 8
        if namelistLineContains(line,"dPhiHatdpsiN")
          line = "dPhiHatdpsiN = " + dPhiHatdpsiNs[i].to_s
        end
        
      end
      
      outFile.write(line + "\n")
    end
    inFile.close
    outFile.close
    
    # Submit job!
    puts "Submitting job #{dirNum}"
    # Make sure job 0 gets submitted first:
    #if i==0
    if true
      # In this way of submitting a job, ruby waits for qsub to complete before moving on.
      puts `cd #{dirNum}; qsub #{jobFilename} &`
    else
      # In this way of submitting a job, ruby does not wait for qsub to complete before moving on.
      job1 = fork do
        exec "cd #{dirNum}; qsub #{jobFilename}"
      end
      Process.detach(job1)
    end
    puts "Done submitting job #{dirNum}"
  end
  
end

puts "Finished submitting jobs for scan."
