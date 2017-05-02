require 'json'
# GeoDiver NameSpace
module GeoDiver
  # module to run the Load the GEO dataset.
  module LoadGeoData
    # To signal error in query sequence or options.
    #
    # ArgumentError is raised when ... exit status is 1; see [1].
    class ArgumentError < ArgumentError
    end

    # To signal internal errors.
    #
    # RuntimeError is raised when there is a problem in writing the input file,
    # running R Script, writing the output etc. These are rare, infrastructure
    # errors, used internally, and of concern only to the admins/developers.
    # One example of a RuntimeError would be R libraries not installed.
    class RuntimeError < RuntimeError
    end

    class << self
      extend Forwardable

      def_delegators GeoDiver, :logger, :public_dir, :db_dir

      # Check if the GEO database has already been downloaded, if not, then
      # download the GEO dataset and extract the meta data and convert into
      # RData
      def run(params, soft_link = true)
        init(params)
        geo_accession  = params['geo_db'].upcase
        meta_json_file = File.join(db_dir, geo_accession,
                                   "#{geo_accession}.json")
        if File.exist? meta_json_file
          logger.debug("Found GeoDb at: '#{meta_json_file}'")
          logger.debug("Parsing GeoDb '#{geo_accession}'")
          meta_data = parse_meta_data(meta_json_file)
        else
          logger.debug("Local GeoDb for '#{geo_accession}' not found.")
          meta_data = download_and_parse_meta_data(geo_accession)
          write_to_json(meta_data, meta_json_file)
        end
        if soft_link
          soft_link_meta_json_to_public_dir(geo_accession, meta_json_file)
        end
        logger.debug('GeoDb loaded into memory')
        meta_data
      end

      def convert_geodb_into_rdata(geo_accession)
        geo_accession = geo_accession.upcase
        return if File.exist?(File.join(db_dir, geo_accession,
                                        "#{geo_accession}.RData"))
        logger.debug("Running: #{load_geo_db_cmd(geo_accession)}")
        Thread.new { system(load_geo_db_cmd(geo_accession)) }
        # TODO: check exit status of the system call
      end

      private

      # Verify paramaters
      def init(params)
        assert_geo_db_present(params)
      end

      #
      def assert_geo_db_present(params)
        logger.debug('Checking if the GEO DB parameter is present.')
        return unless params['geo_db'].nil? || params['geo_db'].empty?
        raise ArgumentError, 'No GEO database provided.'
      end

      def parse_meta_data(meta_json_file)
        logger.debug("Parse the Meta JSON file at: #{meta_json_file}")
        meta_file_content = IO.read meta_json_file
        JSON.parse(meta_file_content)
      end

      #
      def download_and_parse_meta_data(geo_accession)
        file = download_geo_file(geo_accession)
        data = read_geo_file(file)
        return parse_gds_db(data) if geo_accession =~ /^GDS/
        return parse_gse_db(data) if geo_accession =~ /^GSE/
      rescue
        raise ArgumentError, 'GeoDiver was unable to download the GEO Database'
      end

      #
      def download_geo_file(geo_accession)
        if geo_accession =~ /^GDS/
          remote_dir = generate_remote_gds_url(geo_accession)
          file = "#{geo_accession}.soft.gz"
        elsif geo_accession =~ /^GSE/
          remote_dir = generate_remote_gse_url(geo_accession)
          file = "#{geo_accession}_series_matrix.txt.gz"
        end
        output_dir = File.join(db_dir, geo_accession)
        FileUtils.mkdir(output_dir) unless Dir.exist? output_dir
        compressed = File.join(output_dir, file)
        wget_geo_file(remote_dir, compressed, geo_accession, output_dir)
        compressing_geo_file(compressed)
      end

      def wget_geo_file(remote_dir, compressed, geo_accession, output_dir)
        logger.debug("Downloading from: #{remote_dir} ==> #{compressed}")
        `wget #{remote_dir} -O #{compressed} || rm -r #{output_dir}`
        return if $CHILD_STATUS.exitstatus.zero?
        logger.debug "Cannot find Geo Dataset on GEO: #{geo_accession}"
        raise ArgumentError, "Cannot find Geo Dataset on GEO: #{geo_accession}"
      end

      def compressing_geo_file(compressed)
        logger.debug("Uncompressing file: #{compressed.gsub('.gz', '')}")
        system "gunzip --force -c #{compressed} > #{compressed.gsub('.gz', '')}"
        compressed.gsub('.gz', '')
      end

      #
      def generate_remote_gds_url(geo_accession)
        cmd = "bionode-ncbi search gds #{geo_accession} |"\
              " jq -cr 'select(.accession == \"#{geo_accession}\") | .ftplink'"
        url = `#{cmd}`.chomp!
        return if url.nil? || url.empty?
        url + geo_accession + '.soft.gz' 
      end

      def generate_remote_gse_url(geo_accession)
        if geo_accession.length <= 6
          dir = 'GSEnnn'
        else
          if geo_accession.length == 8
            dir_number = geo_accession.match(/GSE(\d\d)\d+/)[1]
          elsif geo_accession.length == 7
            dir_number = geo_accession.match(/GSE(\d)\d+/)[1]
          end
          dir = "GSE#{dir_number}nnn"
        end
        "ftp://ftp.ncbi.nlm.nih.gov/geo/series/#{dir}/#{geo_accession}" \
        "/matrix/#{geo_accession}_series_matrix.txt.gz"
      end

      # Loads the file into memory line by line
      # Stop loading the file once it has read all the meta data.
      def read_geo_file(file)
        data = []
        IO.foreach(file) do |line|
          break if line =~ /^#ID_REF/
          data << line
        end
        data.join
      end

      #
      def parse_gds_db(d)
        {
          'Accession' => d.match(/\^DATASET = (.*)/)[1],
          'Title' => d.match(/!dataset_title = (.*)/)[1],
          'Description' => d.match(/!dataset_description = (.*)/)[1],
          'Sample_Organism' => d.match(/!dataset_platform_organism = (.*)/)[1],
          'Factors' => parse_gds_factors(d),
          'Reference' => d.match(/!Database_ref = (.*)/)[1],
          'Update_Date' => d.match(/!dataset_update_date = (.*)/)[1]
        }
      end

      def parse_gse_db(d)
        {
          'Accession' => d.match(/!Series_geo_accession\t"(.*)"/)[1],
          'Title' => d.match(/!Series_title\t"(.*)"/)[1],
          'Description' => d.match(/!Series_summary\t"(.*)"/)[1],
          'Sample_Organism' => parse_sample_organism(d),
          'Factors' => parse_gse_factors(d),
          'Reference' => d.match(/!Series_relation\t"(.*)"/)[1],
          'Update_Date' => d.match(/!Series_last_update_date\t"(.*)"/)[1]
        }
      end

      #
      def parse_gds_factors(data)
        subsets = data.gsub(/\^DATA.*\n/, '').gsub(/\![dD]ata.*\n/, '')
        factors = {}
        subsets.lines.each_slice(5) do |subset|
          desc = subset[2].match(/\!subset_description = (.*)/)[1]
          type = subset[4].match(/\!subset_type = (.*)/)[1].tr(' ', '.')
          factors[type] ||= {}
          factors[type]['options'] ||= []
          factors[type]['options'] << desc
          factors[type]['value'] = type
        end
        factors
      end

      def parse_gse_factors(data)
        subsets = data.scan(/!Sample_characteristics_ch1\t(.*)/)
        factors = {}
        subsets.each_with_index do |feature, idx|
          a = feature[0].split(/\"?\t?\"/)
          a.delete_if { |e| e =~ /^\s+$/ || e.empty? }
          a.each do |e|
            split = e.split(': ')
            type = split[0]
            factors[type] ||= {}
            factors[type]['value'] = 'characteristics_ch1'
            factors[type]['value'] += ".#{idx}" if idx > 0
            factors[type]['options'] ||= []
            factors[type]['options'] << e
          end
        end
        factors.each { |_, e| e['options'].uniq! }
        factors.delete_if { |_, e| e['options'].size == 1 }
        factors
      end

      def parse_sample_organism(data)
        subset = data.match(/!Sample_organism_ch1\t(.*)/)[1]
        organism = subset.split(/\"?\t?\"/)
        organism.shift
        organism.uniq
      end

      #
      def write_to_json(hash, output_json)
        logger.debug("Writing meta data to file: #{output_json}")
        File.open(output_json, 'w') { |f| f.puts hash.to_json }
      end

      #
      def soft_link_meta_json_to_public_dir(geo_accession, meta_json_file)
        public_meta_json = File.join(public_dir, 'GeoDiver/DBs/',
                                     "#{geo_accession}.json")
        logger.debug("Creating a Soft Link from: #{meta_json_file} ==>" \
                     " #{public_meta_json}")
        return if File.exist? public_meta_json
        FileUtils.ln_s(meta_json_file, public_meta_json)
      end

      #
      def load_geo_db_cmd(geo_accession)
        geo_db_dir = File.join(db_dir, geo_accession)
        "Rscript #{File.join(GeoDiver.root, 'RCore/download_GEO.R')}" \
        " --accession #{geo_accession}" \
        " --outrdata  #{File.join(geo_db_dir, "#{geo_accession}.RData")}" \
        " --geodbDir #{geo_db_dir}" \
        " && echo 'Finished creating Rdata file:" \
        " #{File.join(geo_db_dir, "#{geo_accession}.RData")}'"
      end
    end
  end
end
