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
  sh "cleancss --s0 -s --skip-rebase -o './public/assets/css/style-#{GeoDiver::VERSION}.min.css' './public/assets/css/style.css'"
  sh "cleancss --s0 -s --skip-rebase -o './public/assets/css/home-#{GeoDiver::VERSION}.min.css' './public/assets/css/home.css'"
  sh "uglifyjs './public/assets/js/geodiver.js' './public/assets/js/datatable-materialize.js' './public/assets/js/jquery.filedownload.min.js' -m -c -o './public/assets/js/geodiver-#{GeoDiver::VERSION}.min.js'"
end
