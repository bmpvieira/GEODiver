require 'base64'
require 'English'
require 'forwardable'
require 'json'

require 'geodiver/pool'

# GeoDiver NameSpace
module GeoDiver
  # Module to run the GEO analysis
  module GeoAnalysis
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

      def_delegators GeoDiver, :config, :logger, :public_dir, :users_dir,
                     :db_dir

      #
      def run(params, email, url, load_geo_db_thread)
        init(params, email)
        # wait until geo db has been loaded in background thread
        load_geo_db_thread.join unless load_geo_db_thread.nil?
        run_analysis
        # Compress files in the background
        Thread.new { compress_files(@run_dir, @params['geo_db']) }
        results = generate_results_hash(url)
        save_results_to_file(results)
        results
      end

      private

      def generate_results_hash(url)
        results = { geo_db: @params['geo_db'],
                    user: encode_email,
                    uniq_result_id: @uniq_time,
                    results_url: generate_results_url(url),
                    share_url: generate_share_url(url),
                    assets_path: generate_relative_link(url),
                    full_path: @run_dir,
                    meta_data: parse_meta_json,
                    params: @params }
        add_exit_codes_to_the_results_hash(results)
      end

      def add_exit_codes_to_the_results_hash(results)
        results[:overview_exit_code] = @overview_exit_code
        results[:dgea_exit_code] = @dgea_exit_code
        results[:gage_exit_code] = @gage_exit_code
        results
      end

      #
      def init(params, email)
        @params = params
        @email  = email
        assert_params
        setup_run_and_public_dir
        @uniq_time
      end

      #
      def assert_params
        assert_geo_db_present
        @params['geo_db'] = @params['geo_db'].upcase
        # TODO: assert other Params
      end

      #
      def assert_geo_db_present
        logger.debug('Asserting GEO db is present.')
        return unless @params['geo_db'].nil? || @params['geo_db'].empty?
        fail ArgumentError, 'No GEO database provided.'
      end

      #
      def setup_run_and_public_dir
        @uniq_time = Time.new.strftime('%Y-%m-%d_%H-%M-%S_%L-%N').to_s
        @run_dir = File.join(users_dir, @email, @params['geo_db'], @uniq_time)
        logger.debug("Creating Run Directory: #{@run_dir}")
        FileUtils.mkdir_p(@run_dir)
      end

      def save_results_to_file(results)
        logger.debug("Results: #{results}")
        output_params_json = File.join(@run_dir, 'params.json')
        File.open(output_params_json, 'w') { |f| f.puts results.to_json }
      end

      def run_analysis
        if config[:num_threads] > 1
          run_analysis_multi_threaded
        else
          run_overview
          run_dgea if @params['dgea'] == 'on'
          run_gage if @params['gsea'] == 'on'
        end
      end

      def run_analysis_multi_threaded
        p = Pool.new(config[:num_threads])
        p.schedule { run_overview }
        p.schedule { run_dgea } if @params['dgea'] == 'on'
        p.schedule { run_gage } if @params['gsea'] == 'on'
      ensure
        p.shutdown
      end

      def run_overview
        logger.debug("Running CMD: #{overview_cmd}")
        system(overview_cmd)
        @overview_exit_code = $CHILD_STATUS.exitstatus
      end

      #
      def run_dgea
        logger.debug("Running CMD: #{dgea_cmd}")
        system(dgea_cmd)
        @dgea_exit_code = $CHILD_STATUS.exitstatus
      end

      def run_gage
        logger.debug("Running CMD: #{gsea_cmd}")
        system(gsea_cmd)
        @gage_exit_code = $CHILD_STATUS.exitstatus
      end

      def overview_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/overview.R')}" \
        " --dbrdata #{dbrdata} --rundir '#{@run_dir}/'" \
        " --analyse 'Boxplot,PCA'" \
        " --factor \"#{@params['factor']}\"" \
        " --popA \"#{to_comma_delimited_string(@params['groupa'])}\"" \
        " --popB \"#{to_comma_delimited_string(@params['groupb'])}\"" \
        " --popname1 'Group1' --popname2 'Group2'" \
        ' --dev'
      end

      #
      def dgea_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/dgea.R')}" \
        " --dbrdata #{dbrdata} --rundir '#{@run_dir}/'" \
        " --analyse '#{analyses_to_carry_out.join(',')}'" \
        " --factor \"#{@params['factor']}\"" \
        " --popA \"#{to_comma_delimited_string(@params['groupa'])}\"" \
        " --popB \"#{to_comma_delimited_string(@params['groupb'])}\"" \
        " --popname1 'Group1' --popname2 'Group2'" \
        " --topgenecount #{@params['dgea_number_top_genes']} " \
        ' --foldchange 0 --thresholdvalue 0' \
        " --distance '#{@params['dgea_heatmap_distance_method']}'" \
        " --clustering '#{@params['dgea_heatmap_clustering_method']}'" \
        " --clusterby '#{dgea_clusterby_method}'" \
        " --heatmaprows #{@params['dgea_heatmap_rows']} " \
        " --adjmethod '#{@params['dgea_volcano_pValue_cutoff']}'" \
        " #{option_to_flag(@params['gsea_cluster_by_genes'], '--dendrow')}" \
        " #{option_to_flag(@params['gsea_cluster_by_samples'], '--dendcol')}" \
        ' --dev'
      end

      def gsea_cmd
        "Rscript #{File.join(GeoDiver.root, 'RCore/gage.R')}" \
        " --dbrdata #{dbrdata} --rundir '#{@run_dir}/'" \
        " --factor \"#{@params['factor']}\"" \
        " --popA \"#{to_comma_delimited_string(@params['groupa'])}\"" \
        " --popB \"#{to_comma_delimited_string(@params['groupb'])}\"" \
        " --comparisontype '#{@params['gsea_type']}'"\
        " --genesettype '#{@params['gsea_dataset']}'" \
        " --distance '#{@params['gsea_heatmap_distance_method']}'" \
        " --clustering '#{@params['gsea_heatmap_clustering_method']}'" \
        " --clusterby '#{gage_clusterby_method}'" \
        " --heatmaprows #{@params['gsea_heatmap_rows']} " \
        " #{option_to_flag(@params['gsea_cluster_by_genes'], '--dendrow')}" \
        " #{option_to_flag(@params['gsea_cluster_by_samples'], '--dendcol')}" \
        ' --dev'
      end

      def compress_files(run_dir, geodb)
        cmd = "zip -jr '#{run_dir}/#{geodb}_geodiver_results.zip' '#{run_dir}'"
        logger.debug("Running CMD: #{cmd}")
        system("#{cmd}")
      end

      def analyses_to_carry_out
        analyses = []
        analyses << 'Toptable' if @params['dgea_toptable'] == 'on'
        if @params['dgea_heatmap'] == 'on' || @params['gsea_heatmap'] == 'on'
          analyses << 'Heatmap'
        end
        analyses << 'Volcano' if @params['dgea_volcano'] == 'on'
        analyses
      end

      def to_comma_delimited_string(arr)
        arr.each do |e|
          e.gsub!(/(?<!\\),/, '\,')
          e.gsub!('-', '\-')
        end
        arr.join(',')
      end

      def dgea_clusterby_method
        (@params['dgea_cluster_based_on'] == 'on') ? 'Complete' : 'Toptable'
      end

      def gage_clusterby_method
        (@params['gsea_cluster_based_on'] == 'on') ? 'Complete' : 'Toptable'
      end

      def option_to_flag(option, return_flag)
        return_flag if option == 'on'
      end

      def dbrdata
        File.join(db_dir, @params['geo_db'], "#{@params['geo_db']}.RData")
      end

      def generate_relative_link(url)
        "#{url}/GeoDiver/Users/#{@email}/#{@params['geo_db']}/#{@uniq_time}"
      end

      def generate_results_url(url)
        "#{url}/result/#{encode_email}/#{@params['geo_db']}/#{@uniq_time}"
      end

      def generate_share_url(url)
        "#{url}/sh/#{encode_email}/#{@params['geo_db']}/#{@uniq_time}"
      end

      def parse_meta_json
        meta_json_file = File.join(db_dir, @params['geo_db'],
                                   "#{@params['geo_db']}.json")
        meta_file_content = IO.read meta_json_file
        JSON.parse(meta_file_content)
      end

      def encode_email
        Base64.encode64(@email).chomp
      end
    end
  end
end
