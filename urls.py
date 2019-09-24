from django.urls import path

from . import views

urlpatterns = [
    path('index_of_all_datasets', views.index_all_datasets),
    path('GSE<int:GSE_number>', views.explore_dataset_regex),
    path('GSE<int:GSE_number>_<slug:cell_type>', views.explore_dataset_subtypes),

    path('PRIVATE<int:GSE_number>', views.explore_private_dataset_regex),
    path('PRIVATE<int:GSE_number>_<slug:cell_type>', views.explore_dataset_subtypes),

    path("search_gene_across_all_datasets", views.search_gene_across_all_datasets),
    path("query_gene_expression_all_cells", views.query_gene_expression_all_cells_new),
    # path("query_gene_expression_all_cells_new", views.query_gene_expression_all_cells_new),

    path('usr/cart_ajax', views.submit_order),

    path('scibet_choose_reference', views.scibet_choose_reference),
    path('scibet', views.scibet_general),

    path('00000', views.explore_dataset),

    path('usr_upload', views.explore_usr_upload_dataset),

    path('gene_hint', views.gene_hint),
    path('gene_hint_all_datasets', views.gene_hint_all_datasets),

]