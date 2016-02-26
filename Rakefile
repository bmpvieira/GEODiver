require 'bundler/gem_tasks'
require 'rspec/core'
require 'rspec/core/rake_task'

task default: [:build]

desc 'Builds and installs'
task install: [:installRdependencies, :build] do
  require_relative 'lib/geodiver/version'
  sh "gem install ./geodiver-#{GeoDiver::VERSION}.gem"
end

desc 'Runs tests and builds gem (default)'
task build: [:test] do
  sh 'gem build geodiver.gemspec'
end

desc 'Install R dependencies'
task :installRdependencies do
  sh 'Rscript RCore/Installations.R'
end

task test: :spec
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

task :assets do
  require_relative 'lib/geodiver/version'
  `rm ./public/assets/css/style-*.min.css`
  `rm ./public/assets/css/home-*.min.css`
  `rm ./public/assets/css/app-*.min.css`
  `rm ./public/assets/js/load_css-*.min.js`
  `rm ./public/assets/js/geodiver-*.min.js`
  sh 'cleancss --s0 -s --skip-rebase -o' \
     " './public/assets/css/style-#{GeoDiver::VERSION}.min.css'" \
     " './public/assets/css/style.css'"
  sh 'cleancss --s0 -s --skip-rebase -o' \
     " './public/assets/css/home-#{GeoDiver::VERSION}.min.css'" \
     " './public/assets/css/home.css'"
  sh 'cleancss --s0 -s --skip-rebase -o' \
     " './public/assets/css/app-#{GeoDiver::VERSION}.min.css'" \
     " './public/assets/css/app.css'"
  sh "uglifyjs './public/assets/js/datatable-materialize.js'" \
     " './public/assets/js/jquery.filedownload.min.js'" \
     " './public/assets/js/geodiver.js' -m -c -o" \
     " './public/assets/js/geodiver-#{GeoDiver::VERSION}.min.js'"
  sh "uglifyjs './public/assets/js/load_css.js' -m -c -o" \
     " './public/assets/js/load_css-#{GeoDiver::VERSION}.min.js'"
end
