if (ARGV.length != 2) then
  print "Usage: ruby mkpair.rb <species 1> <species 2>\n"
  exit
end

species1 = ARGV[0]
species2 = ARGV[1]

lines = File.open("primates.config", "r").readlines
lines.each do |line|
  line.chomp!
  data = line.split
  print "<h3>Chromosome #{data[1]}</h3>\n"  
  kv = {}
  dokey = true
  key = ""
  label = ""
  data.each do |dval|
    if (dokey) then
      key = dval
      dokey = false
    else
      kv[key] = dval
      dokey = true
    end
  end
  if (kv[species1] != kv[species2]) then
    label = kv[species2]
  end
  print "<a href=\"./images/#{species1}#{kv[species1]}_#{species2}#{kv[species2]}.txt.png\"><img src=\"./images/#{species1}#{kv[species1]}_#{species2}#{kv[species2]}.x10.png\" /></a>#{label}\n"
end
