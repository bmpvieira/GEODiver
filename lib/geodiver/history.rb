require 'forwardable'

module GeoDiver
  # A class that genrates the history hash
  class History
    extend Forwardable

    def_delegators GeoDiver, :logger, :public_dir, :users_dir, :db_dir

    class <<self
      def determine_history

      end

      def determine_file_structure
        history = {}
        data.each do |d|
          accession, time, file = d
          history[accession] ||= {}
          history[accession][time] ||= []
          history[accession][time] << file
        end
        history
      end

      def all_user_files(user)
        dir = File.join(users_dir, user.info['email'])
        files = Dir.glob(File.join(dir, '**', '*'), File::FNM_DOTMATCH)
        results = []
        files.each do |f|
          next if f =~ /\.$/
          f.gsub!(%r{#{user_dir}/}, '')
          results << f.split('/') unless f.split('/').length < 3
        end
        results
      end
    end
  end
end
