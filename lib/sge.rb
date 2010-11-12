#!/usr/bin/env ruby
#
# = Bio::SGE -- Sun Grid Engine array job submitter (Bio::FlatFile query to SGE)
#
# Copyright::	Copyright (C) 2009, 2010 Toshiaki Katayama <mailto:ktym at hgc dot jp>
# License::	Distributes under the same terms as Ruby
# Site::	http://kanehisa.hgc.jp/~k/sge/
# Download::	http://kanehisa.hgc.jp/~k/sge/sge.rb
# Version::	2.4
#
# == USAGE (AS A COMMAND)
#
# As of the version 2.0, this library can also be used as a command.
#
# Usage:
#     % sge.rb \[options...\] -q input_file -t db_file -c 'command --opts #{query} #{target}'
# 
# Options:
#     -q or --query file
#        Specify a flatfile including multiple entries.
#     -t or --target file
#        Specify a database file to be used.
#     -c or --command 'string'
#        Specify a command line to be executed.
#        The following identifiers can be used in the command line 'string'.
#          '#{query}'       fragmented query file name (== input_file)
#          '#{target}'      target database file name
#          '#{work_dir}'    current working directory
#          '#{task_id}'     SGE_TASK_ID
#          '#{slice}'       -- task_id / @@slice (integer >= 1)
#          '#{input_file}'  -- 'input/#{slice}/#{task_id}'
#          '#{output_file}' -- 'output/#{slice}/#{task_id}'
#          '#{error_file}'  -- 'error/#{slice}/#{task_id}'
#     -o or --sge_opts 'string'
#        Additional options for the qsub command.
#          '-l s_vmem=16G -l mem_req=16' to reserve 16GB RAM for each job
#          '-l cpu_arch=xeon'            to limit to use xeon CPUs only
#        Resource reservation and backfill options:
#          '-R y -l s_rt=12:0:0'         to limit max exec time to 12h (SIGUSER1)
#          '-R y -l h_rt=12:0:0'         to limit max exec time to 12h (SIGKILL)
#          '-R y -pe mpi-fillup 4'       to reserve 4 threads for MPI
#     -m or --task_min integer
#        Start number of tasks (default is 1, increase to start from halfway).
#     -M or --taks_max integer
#        Last value (default is a total number of entries in query).
#     -s or --task_step integer
#        Number of processes per one job (default is 1000). Large value is
#        recommended for short tasks with a large number of queries, and
#        a small value (minimum is 1) can be used for time consuming tasks
#        with a small number of queries.
#     --clear
#        Remove a SGE script and output/error/log directories
#     --clean
#        Remove a count file and the extracted input directory
#     --distclean
#        Exec both of --clear and --clean
#     -h or --help
#        Print this help message.
# 
# Examples:
#     % sge.rb -q data/query.pep -t data/target.pep -c 'blastall -p blastp -i #{query} -d #{target}' -o '-l cpu_arch=xeon'
#     % sge.rb -q data/query.nuc -t /usr/local/db/blast/ncbi/nr -c 'blastall -p blastx -s 10 -i #{query} -d #{target}' -o '-l cpu_arch=xeon -l sjob -l s_vmem=4G,mem_req=4'
#     % sge.rb -q data/dme.nuc -t data/dme.genome -s 1 -c 'exonerate --bestn 1 --model est2genome --showtargetgff 1 --showvulgar yes #{query} #{target}'
#     % sge.rb -q data/hsa.pep -t data/Pfam-A.hmm -m 1000 -M 2000 -s 10 -c 'hmmscan --tblout output/#{slice}/#{task_id}.tbl #{target} #{query}'
#     % sge.rb -q data/refseq.gb -c 'bp_genbank2gff3.pl -out stdout #{query}'
#     % sge.rb --distclean
# 
# See also:
#     http://kanehisa.hgc.jp/~k/sge/
#
# == USAGE (AS A LIBRARY)
# 
# The Bio::SGE class extract entries in a biological flatfile as queries
# and execute a bulk submission to the Sun Grid Engine as an array job.
#
# This class takes a flatfile (e.g. multi FASTA file) as a 'query',
# a database file as a 'target', and a command line to be executed
# as a 'command' (see also SCRIPT VARIABLES section).
#
# The flatfile must be accepted by the Bio::FlatFile.auto class method
# of the BioRuby (http://bioruby.org/) package.
#
# Instantiation of the Bio::SGE object can be done by
#
#   sge = Bio::SGE.new(query, target, command, sge_opts)
# 
# or by assigning these values through accessors prior to a job submission
#
#   sge = Bio::SGE.new
#   sge.query = 'flat_file'
#   sge.target = 'target_database_file'
#   sge.command = 'command --to_be_executed --with_opts'
#
# or by assigning these values with a block parameter.
#
#   sge = Bio::SGE.new { |opt|
#     opt.query = 'flat_file'
#     opt.target = 'target_database_file'
#     opt.command = 'command --to_be_executed --with_opts'
#   }
#
# Then, the "prepare" method will
# 
# * create output directories
# * generate a SGE script to be submitted
# * extract each entry in the query as separate files
#   (files are numbered by the order of appearance)
# 
# and now you can submit your SGE job by the "submit" method.
# 
#   sge.prepare
#   sge.submit
# 
# The "submit" method will automatically take care of messy tasks such that
# (1) splitting array jobs according to the number of total jobs, (2) save
# stdout and stderr from SGE system to a separate log directory etc.
# 
# == RESULTS
# 
# The execution results will be stored in the following files and directories.
# 
#   count.txt     # correspondence table of the file numbers and entry IDs
#   input/        # extracted sequence files (one file, one sequence)
#   output/       # outputs of the command (numberd same as the input files)
#   error/        # errors of the command (numberd same as the input files)
#   log/          # log files of the qsub run (stdout and stderr)
# 
# You can confirm whether there were no system errors during the SGE execution
# by sizes and contents of files in the log/ directory.
# 
# Then, check the error/ directory whether there was a problem or not in your
# jobs (some command may utilize the stderr to another purpose).
# 
# Finally, main results can be obtained from files in the output/ directory.
# 
# == ADVANCED USAGE
# 
# You can individually call following methods instead of the "prepare" method.
# 
#   sge.setup     # to prepare output directories
#   sge.script    # to generate a SGE script
#   sge.extract   # to extract each entry
# 
# Therefore, if you want to reuse the sequence files already extracted to
# the input directory, just comment out the line calling "prepare" method
# (and also avoid to use "extract" method, of course).
# 
#   #sge.prepare  # comment out this line in your script
#   sge.script
#   sge.setup
#   #sge.extract  # don't use this as well
# 
#   sge.submit    # then submit
# 
# Reversely, you can also clean up the working directory (e.g. to remove
# test or previous execution results) by the following methods.
# 
#   sge.clear     # to remove a SGE script and output/error/log directories
#   sge.clean     # to remove a count file and the extracted input directory
#   sge.distclean # to remove all of the above
# 
# == SGE OPTIONS
# 
# You can specify the "-t start-last:step" range values for a array job
# by following accessors (these are optional; see EXAMPLES section below).
# 
#   sge.task_min  # start value (default is 1)
#   sge.task_max  # last value (default is a total number of entries in query)
#   sge.task_step # number of processes per one job (default is 1000)
#   sge.sge_opts  # additional options for the qsub command
# 
# For example, if you only need to calculate on sequences starting from 8421st
# upto 9064th, and want to invoke 100 processes per each qsub execution, you
# can specify them by the following way.
# 
#   sge.task_min = 8421
#   sge.task_max = 9064
#   sge.task_step = 100
#   
#   sge.submit
# 
# == OVER ALL SKELETON
# 
#   #!/usr/bin/env ruby
# 
#   require 'sge'
#
#   sge = Bio::SGE.new { |opt|
#     opt.query = 'flat_file'
#     opt.target = 'target_database_file'
#     opt.command = 'command --to_be_executed --with_opts'
#     opt.sge_opts = '-l cpu_arch=xeon'
#     opt.task_min = 8421
#     opt.task_max = 9064
#     opt.task_step = 100
#   }
#   sge.clear           # included in sge.distclean
#   sge.clean		# included in sge.distclean
#   sge.script          # included in sge.prepare
#   sge.setup           # included in sge.prepare
#   sge.extract         # included in sge.prepare
#   sge.submit
# 
# == SCRIPT VARIABLES
# 
# In the 'command' specification, you can use following identifiers as variables.
# 
#   '#{query}'       fragmented query file name (== input_file)
#   '#{target}'      target database file name
#   '#{work_dir}'    current working directory
# 
#   '#{task_id}'     SGE_TASK_ID
#   '#{slice}'       -- task_id / @@slice (integer >= 1)
#   '#{input_file}'  -- 'input/#{slice}/#{task_id}'
#   '#{output_file}' -- 'output/#{slice}/#{task_id}'
#   '#{error_file}'  -- 'error/#{slice}/#{task_id}'
# 
# Note that these identifires must be kept in 'single quotes' to avoid variable
# expansion before the script generation (see EXAMPLES section in below).
# 
# == EXAMPLES
# 
# 1. Exonerate (Query = Multi fasta protein sequences; Target = Genomic DNA)
# 
#     #!/usr/bin/env ruby
# 
#     require 'sge'
# 
#     sge = Bio::SGE.new { |opt|
#       opt.query = 'd.melanogaster.pep'
#       opt.target = 'genomic_scaffolds'
#       opt.command = 'exonerate --bestn 1 --model protein2genome --showtargetgff 1 --showvulgar yes #{query} #{target}'
#       opt.sge_opts = '-l cpu_arch=xeon'
#     }
#     sge.prepare
#     sge.submit
# 
# 
# 2. BLAST (Query = Multi fasta; Target = BLAST DB)
# 
#     #!/usr/bin/env ruby
# 
#     require 'sge'
# 
#     sge = Bio::SGE.new { |opt|
#       opt.query = 'query.pep'
#       opt.target = 'target.pep'
#       opt.command = 'blastall -p blastp -i #{query} -d #{target}'
#       opt.sge_opts = '-l cpu_arch=xeon'
#     }
#     sge.prepare
#     sge.submit
# 
# 
# 3. HMMER (Query = Multi fasta protein sequences; Target = Pfam DB)
# 
#     #!/usr/bin/env ruby
# 
#     require 'sge'
# 
#     sge = Bio::SGE.new { |opt|
#       opt.query = 'data/h.sapiens.pep'
#       opt.target = 'db/Pfam_ls'
#       opt.command = 'hmmscan --tblout output/#{slice}/#{task_id}.tbl #{target} #{query}'
#       opt.sge_opts = '-l cpu_arch=xeon'
#     }
#     sge.prepare
#     sge.submit
#
# 4. RefSeq to GFF (Query = RefSeq entries in GenBank format)
#
#     #!/usr/bin/env ruby
# 
#     require 'sge'
# 
#     sge = Bio::SGE.new { |opt|
#       opt.query = 'invertebrate6.genomic.gbff'
#       opt.command = 'bp_genbank2gff3.pl -out stdout #{query}'
#     }
#     sge.prepare
#     sge.submit
#
# == CHANGE LOG
#
# === 2010/09/11 v2.4
#
# * changed to skip "extract" when "count.txt" file exists, so that
#   user can easily re-submit the job (with different parameter or fix)
#   just after the --clear. To extract again (starting from scratch),
#   use --clean (with --clear) or --distclean first.
# * doc fix
#
# === 2010/05/21 v2.3
#
# * slice is changed to start from 1 (instead of 0) and have 1000 files
#   per directory (instead of 10000).
#
# === 2010/03/25 v2.2
#
# * doc fix
#
# === 2009/12/10 v2.1
#
# * --clear, --clean, --distcrean options are supported.
#
# === 2009/12/07 v2.0
#
# * Extended to be used as a command.
#
# === 2009/11/13 v1.3
#
# * SGE class is moved under the Bio name space (Bio::SGE) as it tightly
#   depends on the Bio::FlatFile (in BioRuby).
# * Bio::SGE is improved to accept options as a block parameter as well.
#
# === 2009/11/02 v1.2
#
# * slice functionality is fixed to properly create slice directories
#   under the output and error directories
#
# === 2009/09/29 v1.1
# * slice (sub directory to supress the number of files in a single directory)
#   is introduced not to overload file server (MDS)
# * fixed document
# 
# === 2009/09/29 v1.0
# * SGE_TASK_LAST is introduced to avoid remainder jobs are submitted
# * documentation is rewrited in Rdoc format
# * web site is opend and released to public
# 
# === 2009/09/29 v0.3
# * SGE_TASK_STEPSIZE is introduced not to overload the SGE manager
#   by a bunch of short time jobs
# 
# === 2009/09/23 v0.2
# * query/target variables are intoduced to allow commands having
#   BLAST-like options for specifying query and target files
# * added documentation
# 
# === 2009/09/17 v0.1
# * implemented FASTA file extraction and qsub submission functionality
#

require 'rubygems'
require 'fileutils'
require 'bio'

module Bio

class SGE

  # Number of files per directory
  @@slice = 1000

  # Template string for script generation
  @@template = <<'END'
#$ -S /usr/local/bin/ruby

work_dir = "%WORK_DIR%"

offset = ENV["SGE_TASK_ID"].to_i
limit  = ENV["SGE_TASK_STEPSIZE"].to_i
last   = ENV["SGE_TASK_LAST"].to_i

slice = slice_old = nil

offset.upto(offset + limit) do |task_id|
  break if task_id > last

  slice_old = slice
  slice = (task_id - 1) / %SLICE% + 1
  output_dir = "%OUTPUT_DIR%/#{slice}"
  error_dir = "%ERROR_DIR%/#{slice}"
  Dir.mkdir(output_dir) if slice_old != slice and ! File.directory?(output_dir)
  Dir.mkdir(error_dir)  if slice_old != slice and ! File.directory?(error_dir)

  input_file  = "%INPUT_DIR%/#{slice}/#{task_id}"
  output_file = "%OUTPUT_DIR%/#{slice}/#{task_id}"
  error_file  = "%ERROR_DIR%/#{slice}/#{task_id}"

  query = input_file
  target = "%TARGET%"

  if File.exists?(query)
    system("%COMMAND% > #{output_file} 2> #{error_file}")
  end
end
END

  attr_accessor :query, :target, :command, :sge_opts, :count
  attr_accessor :task_min, :task_max, :task_step
  attr_accessor :work_dir, :log_dir, :input_dir, :output_dir, :error_dir

  def initialize(query = nil, target = nil, command = nil, sge_opts = nil)
    @work_dir = Dir.pwd
    @query = "#{@work_dir}/#{query}"
    @target = "#{@work_dir}/#{target}"
    @command = command
    @sge_opts = sge_opts

    yield(self) if block_given?

    @log_dir = "log"
    @input_dir = "input"
    @output_dir = "output"
    @error_dir = "error"
    @script_file = "script.rb"
    @count_file = "count.txt"
  end

  def prepare
    setup
    script
    extract
  end

  def submit
    unless @count
      $stderr.puts "Reading #{@count_file} ..."
      @count = File.readlines(@count_file).last[/^\d+/].to_i
      $stderr.puts "done."
    end

    task_min = @task_min || 1
    task_max = @task_max || @count
    task_step = @task_step || 1000

    # system upper limit is 75000
    limit = 50000
    task_min.step(task_max, limit) do |offset|
      opts = "#{@sge_opts} -o #{@log_dir} -e #{@log_dir} -cwd"
      span = "-t #{offset}-#{[offset + limit, task_max].min}:#{task_step}"
      qsub = "qsub #{opts} #{span} #{@script_file}"
      $stderr.puts "Submitting ... #{qsub}"
      system(qsub)
    end
  end

  def rmtree(file)
    $stderr.print "Deleting #{file} ... "
    FileUtils.rmtree(file)
    $stderr.puts "done."
  end

  def clear
    rmtree(@script_file)
    rmtree(@output_dir)
    rmtree(@error_dir)
    rmtree(@log_dir)
  end

  def clean
    rmtree(@count_file)
    rmtree(@input_dir)
  end

  def distclean
    clear
    clean
  end

  def mkpath(dir)
    $stderr.print "Creating #{dir} ... "
    if File.directory?(dir)
      $stderr.puts "skip (already exists)."
    else
      FileUtils.mkpath(dir)
      $stderr.puts "done."
    end
  end

  def setup
    mkpath(@log_dir)
    mkpath(@input_dir)
    mkpath(@output_dir)
    mkpath(@error_dir)
  end

  def script
    sge_script = @@template.dup
    sge_script.gsub!('%WORK_DIR%', @work_dir)
    sge_script.gsub!('%INPUT_DIR%', @input_dir)
    sge_script.gsub!('%OUTPUT_DIR%', @output_dir)
    sge_script.gsub!('%ERROR_DIR%', @error_dir)
    sge_script.gsub!('%TARGET%', @target)
    sge_script.gsub!('%COMMAND%', @command)
    sge_script.gsub!('%SLICE%', @@slice.to_s)

    File.open(@script_file, "w") do |file|
      file.puts sge_script
    end
  end

  def extract
    return if File.exists?(@count_file)

    slice = slice_old = nil
    @count = 0
    File.open(@count_file, "a") do |count_file|
      Bio::FlatFile.auto(@query) do |ff|
        ff.each do |entry|
          @count += 1
          $stderr.print "Extracting ... #{@count} (#{entry.entry_id}) "
          if (@task_min and @count < @task_min) or (@task_max and @count > @task_max)
            $stderr.puts "skip."
            next
          else
            slice_old = slice
            slice = (@count - 1) / @@slice + 1
            slice_dir = "#{@input_dir}/#{slice}"
            mkpath(slice_dir) if slice_old != slice
            File.open("#{slice_dir}/#{@count}", "w") do |file|
              file.puts ff.entry_raw
            end
            count_file.puts [@count, entry.entry_id].join("\t")
            $stderr.puts "done."
          end
        end
      end
    end
  end


end # class SGE

end # module Bio


if __FILE__ == $0

  def show_usage
    prog  = File.basename($0)
    usage = %Q[
Usage:
    % #{prog} \[options...\] -q input_file -t db_file -c 'command --opts \#{query} \#{target}'

Options:
    -q or --query file
       Specify a flatfile including multiple entries.
    -t or --target file
       Specify a database file to be used.
    -c or --command 'string'
       Specify a command line to be executed.
       The following identifiers can be used in the command line 'string'.
         '\#{query}'       fragmented query file name (== input_file)
         '\#{target}'      target database file name
         '\#{work_dir}'    current working directory
         '\#{task_id}'     SGE_TASK_ID
         '\#{slice}'       -- task_id / @@slice (integer >= 0)
         '\#{input_file}'  -- "input/\#{slice}/\#{task_id}"
         '\#{output_file}' -- "output/\#{slice}/\#{task_id}"
         '\#{error_file}'  -- "error/\#{slice}/\#{task_id}"
    -o or --sge_opts 'string'
       Additional options for the qsub command.
         '-l s_vmem=16G -l mem_req=16' to reserve 16GB RAM for each job
         '-l cpu_arch=xeon'            to limit to use xeon CPUs only
       Resource reservation and backfill options:
         '-R y -l s_rt=12:0:0'         to limit max exec time to 12h (SIGUSER1)
         '-R y -l h_rt=12:0:0'         to limit max exec time to 12h (SIGKILL)
         '-R y -pe mpi-fillup 4'       to reserve 4 threads for MPI
    -m or --task_min integer
       Start number of tasks (default is 1, increase to start from halfway).
    -M or --taks_max integer
       Last value (default is a total number of entries in query).
    -s or --task_step integer
       Number of processes per one job (default is 1000). Large value is
       recommended for short tasks with a large number of queries, and
       a small value (minimum is 1) can be used for time consuming tasks
       with a small number of queries.
    -h or --help
       Print this help message.
    --clear
       Remove a SGE script and output/error/log directories
    --clean
       Remove a count file and the extracted input directory
    --distclean
       Exec both of --clear and --clean

Examples:
    % #{prog} -q data/query.pep -t data/target.pep -c 'blastall -p blastp -i \#{query} -d \#{target}' -o '-l cpu_arch=xeon'
    % #{prog} -q data/query.nuc -t /usr/local/db/blast/ncbi/nr -c 'blastall -p blastx -s 10 -i \#{query} -d \#{target}' -o '-l cpu_arch=xeon -l sjob -l s_vmem=4G,mem_req=4'
    % #{prog} -q data/dme.nuc -t data/dme.genome -s 1 -c 'exonerate --bestn 1 --model est2genome --showtargetgff 1 --showvulgar yes \#{query} \#{target}'
    % #{prog} -q data/hsa.pep -t data/Pfam-A.hmm -m 1000 -M 2000 -s 10 -c 'hmmscan --tblout output/\#{slice}/\#{task_id}.tbl \#{target} \#{query}'
    % #{prog} -q data/refseq.gb -c 'bp_genbank2gff3.pl -out stdout \#{query}'
    % #{prog} --distclean

See also:
    http://kanehisa.hgc.jp/~k/sge/

]
    puts usage
    exit
  end

  require 'getoptlong'

  $opts = Hash.new

  args = GetoptLong.new(
    [ '--query',     '-q',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--target',    '-t',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--command',   '-c',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--sge_opts',  '-o',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--task_min',  '-m',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--task_max',  '-M',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--task_step', '-s',  GetoptLong::REQUIRED_ARGUMENT ],
    [ '--clear',            GetoptLong::NO_ARGUMENT ],
    [ '--clean',            GetoptLong::NO_ARGUMENT ],
    [ '--distclean',        GetoptLong::NO_ARGUMENT ],
    [ '--help',      '-h',  GetoptLong::NO_ARGUMENT ]
  )

  args.each_option do |name, value|
    case name
    when /--query/
      $opts[:query] = value
    when /--target/
      $opts[:target] = value
    when /--command/
      $opts[:command] = value
    when /--sge_opts/
      $opts[:sge_opts] = value
    when /--task_min/
      $opts[:task_min] = value.to_i
    when /--task_max/
      $opts[:task_max] = value.to_i
    when /--task_step/
      $opts[:task_step] = value.to_i
    when /--clear/
      $opts[:clear] = true
    when /--clean/
      $opts[:clean] = true
    when /--distclean/
      $opts[:clear] = true
      $opts[:clean] = true
    when /--help/
      $opts[:help] = true
    end
  end

  if $opts[:clear]
    sge = Bio::SGE.new
    sge.clear
  end

  if $opts[:clean]
    sge = Bio::SGE.new
    sge.clean
  end

  show_usage if $opts[:help] or !$opts[:command]

  sge = Bio::SGE.new { |opt|
    opt.query     = $opts[:query]     if $opts[:query]
    opt.target    = $opts[:target]    if $opts[:target]
    opt.command   = $opts[:command]   if $opts[:command]
    opt.sge_opts  = $opts[:sge_opts]  if $opts[:sge_opts]
    opt.task_min  = $opts[:task_min]  if $opts[:task_min]
    opt.task_max  = $opts[:task_max]  if $opts[:task_max]
    opt.task_step = $opts[:task_step] if $opts[:task_step]
  }
  sge.prepare
  sge.submit
end
