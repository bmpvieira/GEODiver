require 'fileutils'
require 'forwardable'
require 'geodiver/load_geo_db'
require 'geodiver/pool'

# GeoDiver NameSpace
module GeoDiver
  # Module to create the Cache
  module GdsCache
    class <<self
      extend Forwardable

      def_delegators GeoDiver, :logger

      def run
        logger.debug "GeoDiver DB Directory #{GeoDiver.db_dir}"
        upperlimit = check_gds_ftp_folders * 1
        p = Pool.new(GeoDiver.config[:num_threads])
        (1..upperlimit).each do |idx|
          p.schedule(idx) { |i| generate_cache("GDS#{i}") }
        end
      ensure
        p.shutdown
      end

      def generate_cache(geo_accession)
        params = { 'geo_db' => geo_accession }
        logger.debug "#{geo_accession}: Downloading and extracting meta data"
        LoadGeoData.run(params, false)
        cmd = load_geo_db_cmd(geo_accession)
        logger.debug "#{geo_accession}: Running: #{cmd}"
        `#{cmd}`
        logger.debug "#{geo_accession}: Exit Code: #{$CHILD_STATUS.exitstatus}"
        if $CHILD_STATUS.exitstatus == 0
          logger.debug "#{geo_accession}: Cleaning up folder"
          cleanup(geo_accession)
        else
          FileUtils.touch File.join(GeoDiver.db_dir, geo_accession, '.failed')
        end
      rescue ArgumentError
        logger.debug "GEO DB not found: #{geo_accession} "
      end

      def load_geo_db_cmd(geo_accession)
        geo_db_dir = File.join(GeoDiver.db_dir, geo_accession)
        log_file = File.join(geo_db_dir, 'log_file.txt')
        cmd = "Rscript #{File.join(GeoDiver.root, 'RCore/download_GEO.R')}" \
              " --accession #{geo_accession}" \
              " --outrdata #{File.join(geo_db_dir, "#{geo_accession}.RData")}" \
              " --geodbDir #{geo_db_dir} &> #{log_file}"
        "echo #{cmd} && #{cmd}"
      end

      def check_gds_ftp_folders
        return 20
        i = 1
        loop do
          `wget -q --spider ftp://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS#{i}nnn/`
          exit_code = $CHILD_STATUS.exitstatus
          break if exit_code != 0
          i += 1
        end
        i
      end

      def cleanup(geo_accession)
        soft_file = File.join(GeoDiver.db_dir, geo_accession, '*.soft')
        `rm #{soft_file}`
        `rm #{soft_file}.gz`
      end
    end
  end
end
