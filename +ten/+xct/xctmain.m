   % # @profile(precision=4)
   %  def get_embeds(self,
   %               train = True;
   %               n_steps = 1000;
   %               lr = 0.01;
   %               plot_losses: bool = False,
   %               losses_file_name: str = None,
   %               dist_metric: str = "euclidean",
   %               rank: bool = False,
   %               **optim_kwargs
   %               ):
   %      if train:
   %          projections = self._nn_trainer.train(n_steps=n_steps, 
   %                                                  lr=lr; 
   %                                                  **optim_kwargs)
   %      else:
   %          projections = self._nn_trainer.reload_embeds() # reload projections
   %      if plot_losses:
   %          self._nn_trainer.plot_losses(losses_file_name)
   % 
   %      self._aligned_result = self._nn_trainer.nn_aligned_dist(projections,
   %                                                              gene_names_x=self._genes[self._cell_names[0]],
   %                                                              gene_names_y=self._genes[self._cell_names[1]],
   %                                                              w12_shape=self.w12_shape;
   %                                                              dist_metric=dist_metric;
   %                                                              rank=rank)
   %      return projections
   % 
   % 
   % 
   % 
   % 
   % 
   % 
   % 
   % 
   %      % ------------------------
   % 
   %              # cal w
   %      self._w, self.w12_shape = self._build_w(alpha=alpha,
   %                                              query_DB=query_DB;
   %                                              scale_w=scale_w;
   %                                              mu=mu)
   % 
   %      self._nn_trainer = ManifoldAlignmentNet(self._get_data_arr(),
   %                                              w=self._w,
   %                                              n_dim=n_dim;
   %                                              layers=None)