function lattice_visualization(visualize_lattice,lattice)
% LATTICE_VISUALIZATION Display the lattice.

if visualize_lattice == 1
    if length(unique(lattice)) == 1
        lattice_visualization = rescale(lattice,1,1);
    else
        lattice_visualization = rescale(lattice,0,1);
    end
    
    lattice_visualization = fliplr(rot90(lattice_visualization,-1));
    imshow(lattice_visualization)
    
    saveas(gcf,[pwd '/temp_results/lattice/lattice.png']);
end

end