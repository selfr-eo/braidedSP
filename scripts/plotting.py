

def plot_smoothing():

    if savePlot == True:
        # Step 5. Find endpoints
        endpoints = find_endpoints(skeleton)

        # display results and save to figure
        if second_dilation == True:
            fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(16, 4), sharex=True, sharey=True)
        else:
            fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(16, 4), sharex=True, sharey=True)
        ax = axes.ravel()

        ax[0].imshow(mask_image, cmap=plt.cm.gray)
        ax[0].axis('off')
        ax[0].set_title('Water mask', fontsize=15)

        ax[1].imshow(dilated_mask, cmap=plt.cm.gray)
        ax[1].axis('off')
        ax[1].set_title('Dilation', fontsize=15)

        ax[2].imshow(largest_component_mask, cmap=plt.cm.gray)
        ax[2].axis('off')
        ax[2].set_title('Largest connection', fontsize=15)

        ax[3].imshow(smoothed_mask, cmap=plt.cm.gray)
        ax[3].axis('off')
        ax[3].set_title('Gaussian filter', fontsize=15)

        if second_dilation == True:
            ax[4].imshow(largest_component_mask_smoothed, cmap=plt.cm.gray)
            ax[4].axis('off')
            ax[4].set_title('Secondary dilation', fontsize=15)

            ax[5].imshow(np.zeros_like(skeleton), cmap='gray')
            y, x = np.where(skeleton == 1)  # Get skeleton coordinates
            ax[5].scatter(x, y, color='red', s=0.05)  # Plot skeleton points as red dots
            ax[5].axis('off')
            ax[5].set_title('Skeleton', fontsize=15)

        else:
            ax[4].imshow(np.zeros_like(skeleton), cmap='gray')
            y, x = np.where(skeleton == 1)  # Get skeleton coordinates
            ax[4].scatter(x, y, color='red', s=0.05)  # Plot skeleton points as red dots
            ax[4].axis('off')
            ax[4].set_title('Skeleton', fontsize=15)

        fig.tight_layout()
        #fig.suptitle("Dilation: "+str(dilate_amount)+", Gauss: "+str(gauss_amount))
        plt.savefig(figdir+'/'+pixcdate+'gauss_dilate_'+str(gauss_amount)+str(dilate_amount)+'extractedSkeleton.png')
        plt.close()


def plot_labeled_skeleton():

    if savePlot == True:
    # Get coordinates and labels of all labeled pixels
    y_coords, x_coords = np.where(labeled_skeleton > 0)
    labels = labeled_skeleton[y_coords, x_coords]

    # display results
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(4, 4), sharex=True, sharey=True)

    ax = axes.ravel()

    ax[0].imshow(np.zeros_like(final_skeleton), cmap='gray')
    y, x = np.where(final_skeleton == 1)  # Get skeleton coordinates
    ax[0].scatter(x, y, color='red', s=0.05)  # Plot skeleton points as red dots
    y, x = np.where(endpoints == 1) 
    ax[0].scatter(x, y, color='yellow', s=0.5) 
    y, x = np.where(joints == 1) 
    ax[0].scatter(x, y, color='blue', s=0.5) 
    ax[0].axis('off')
    ax[0].set_title('E&J', fontsize=15)


    # View label
    ax[1].imshow(np.zeros_like(final_skeleton), cmap='gray')
    ax[1].scatter(x_coords, y_coords, c=labels, cmap='tab20', s=1)

    # Annotate each branch with its ID near its center
    unique_labels = np.unique(labels)
    for branch_id in unique_labels:
        branch_yx = np.column_stack(np.where(labeled_skeleton == branch_id))
        if len(branch_yx) > 0:
            # Find approximate center by using the median position
            center_y, center_x = np.median(branch_yx, axis=0).astype(int)
            ax[1].text(center_x+300, center_y, str(branch_id), color="white", fontsize=8, ha='center')

    ax[1].axis('off')
    y, x = np.where((final_skeleton == 1) & (labeled_skeleton == 0)) 
    ax[1].scatter(x, y, color='yellow', s=1) # These are skeleton points without a label


    ax[1].set_title('Labeled branch', fontsize=15)

    fig.tight_layout()
    plt.savefig(figdir+'/'+pixcdate+'labeledSkeleton.png')
    plt.close()


def plot_merged(cl_merged, maskdate,figdir):

    # Ensure GeoDataFrame is in a projected CRS for accurate placement
    if not cl_merged.crs.is_projected:
        # print("Reprojecting to UTM for better plotting.")
        cl_merged = cl_merged.to_crs('EPSG:3857')  # Example of a projected CRS

    # Create a color map for branch IDs
    unique_ids = cl_merged['branch_id'].unique()
    num_branches = len(unique_ids)
    # Generate a repeating discrete color map using tab20 or a larger palette
    colors = plt.cm.tab20(np.linspace(0, 1, 20))  # tab20 has 20 colors
    color_map = {branch_id: colors[i % 20] for i, branch_id in enumerate(unique_ids)}

    # Plot the linestrings
    fig, ax = plt.subplots(figsize=(3, 6))
    for _, row in cl_merged.iterrows():
        color = color_map[row['branch_id']]
        ax.plot(*row.geometry.xy, color=color, linewidth=2, label=f"Branch {row['branch_id']}")

    # Remove duplicate labels in the legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    #ax.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1, 1))

    # Add text labels for each branch ID
    for _, row in cl_merged.iterrows():
        x, y = row.geometry.xy[0][len(row.geometry.xy[0]) // 2], row.geometry.xy[1][len(row.geometry.xy[1]) // 2]
        ax.text(x, y, str(row['branch_id']), fontsize=5, ha='center', va='center', color='black')

    # Set plot properties
    ax.set_title("Generated "+str(maskdate))
    ax.axis("equal")
    # Remove x and y tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(figdir+'/'+maskdate+'_generated_centerlines.png')
    # plt.show()
    plt.close()