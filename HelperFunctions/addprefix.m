function str = addprefix(p,args)

    str = reshape(cat(1,cellfun(@(s)[p,'.',s],args(1:2:end-1),...
              'UniformOutput',false),args(2:2:end)),1,[]);

end